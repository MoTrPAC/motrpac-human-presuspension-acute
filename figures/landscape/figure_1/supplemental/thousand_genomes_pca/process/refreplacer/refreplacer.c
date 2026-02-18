#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <dirent.h>
#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#define N_ARGS 6
#define ARG_MERGED_VCF_PATH 1
#define ARG_SAMPLE_DIRECTORY_PATH 2
#define ARG_THREADS 3
#define ARG_CHROM 4
#define ARG_OUTPUT_VCF_PATH 5
#define PATH_BUFFER 1024
#define REF_ALLELE 2
#define MISSING_ALLELE 0

/**
 * Initializes the sync reader. The sync reader is set up to pair matching VCF coordinates across
 * all files. If multiple threads are supported, they are configured for the sync reader. The
 * VCF file merged by the BIC is added to the reader, and all VCF files in the sample directory are
 * also added.
 *
 * @param sr The sync reader from htslib.
 * @param merged_vcf_path The path to the merged VCF file from the BIC.
 * @param sample_directory_path The path to the directory containing the sample VCF files.
 * @param threads The number of threads to use.
 * @param chrom The chromosome being parsed.
 */
void initialize_sync_reader(bcf_srs_t *sr, char *merged_vcf_path, char *sample_directory_path, int n_threads, char *chrom) {

    // Automatically pair loci with matching reference alleles
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);

    // Set the number of threads to use
    if (n_threads > 1) {
        bcf_sr_set_threads(sr, n_threads - 1);
    }

    // Set chromosome as current region
    bcf_sr_set_regions(sr, chrom, false);

    // Add files to the sync reader
    bcf_sr_add_reader(sr, merged_vcf_path);

    // Get directory from path
    DIR *dir = opendir(sample_directory_path);
    if (dir == NULL) {
        fprintf(stderr, "Directory containing samples cannot be opened.");
        exit(EXIT_FAILURE);
    }

    // Read each entry in the directory
    struct dirent *dir_entry;
    while ((dir_entry = readdir(dir)) != NULL) {

        // Check if entry ends in .vcf.gz
        if (strcmp(dir_entry->d_name + strlen(dir_entry->d_name) - 7, ".vcf.gz") == 0) {

            // Add VCF sample file to reader
            char path[PATH_BUFFER];
            strcpy(path, sample_directory_path);
            strcat(path, "/");
            strcat(path, dir_entry->d_name);
            bcf_sr_add_reader(sr, path);
        }
    }
    closedir(dir);
}

/**
 * Creates a map from the sample files in the sync reader to the appropriate column in the output
 * VCF file.
 *
 * @param sr The sync reader from htslib.
 * @param sample_indices An array to populate with the indices representing matched columns.
 * @param out_header The header of the output VCF file to match sample names with.
 */
void map_readers_to_samples(bcf_srs_t *sr, int *sample_indices, bcf_hdr_t *out_header) {

    // Iterate over each reader in the sync reader
    for (int n = 0; n < sr->nreaders; n++) {

        // For each reader, iterate over each column name in the output header
        for (int m = 0; m < sr->nreaders; m++) {

            // Check if the output header's sample name matches the current reader's sample name
            if (strcmp(out_header->samples[m], sr->readers[n].header->samples[0]) == 0) {
                sample_indices[n] = m;
            }
        }
    }
}

/**
 * Get the index of GT in the format field of the specified BCF line.
 *
 * @param line The line to check for the GT index.
 * @param fmt The format field of the current line.
 * @param sample_header The header of the associated reader.
 * @return The index of GT in the format field. Will return -1 if not found.
 */
int get_gt_index(bcf1_t *line, bcf_fmt_t *fmt, bcf_hdr_t *sample_header) {

    // Ensure that the format field is present
    if (line->n_sample && line->n_fmt) {

        // Iterate over each format field
        for (int l = 0; l < line->n_fmt; l++) {

            // Check if specification in the header matches GT
            if (strcmp(sample_header->id[BCF_DT_ID][fmt[l].id].key, "GT") == 0) {
                return l;
            }
        }
    }

    return -1;
}

/**
 * Tests if the specified sample on the VCF record is homozygous for a specific allele.
 *
 * @param fmt The format field of the current line.
 * @param gt_idx The index of GT in the format field.
 * @param sample_idx The index of the sample to check the genotype of.
 * @param value The value of the allele to test for.
 * @return True if the sample is homozygous for the allele specified.
 */
bool is_homozygous_for(bcf_fmt_t *fmt, int gt_idx, int sample_idx, int value) {

    // Generate a pointer to the sample's data in the VCF record
    u_int8_t *ptr = (fmt->p + sample_idx * fmt->size);

    // Generalize the notion of conversion and value testing for various integer sizes
    #define IS_HOM_SWITCH_PATH(converter, type) return (converter(ptr) == value) && (converter(ptr + sizeof(type)) == value);

    // Check which integer size is used for the variant, and convert data accordingly
    switch (fmt[gt_idx].type) {
        case BCF_BT_INT8: IS_HOM_SWITCH_PATH(le_to_i8, int8_t); break;
        case BCF_BT_INT16: IS_HOM_SWITCH_PATH(le_to_i16, int16_t); break;
        case BCF_BT_INT32: IS_HOM_SWITCH_PATH(le_to_i32, int32_t); break;
    }
}

/**
 * Write a homozygous genotype of a specified value to the VCF line for a specified sample.
 *
 * @param fmt The format field of the current line.
 * @param gt_idx The index of GT in the format field.
 * @param sample_idx The index of the sample to write the genotype to.
 * @param value The value to write.
 */
void write_homozygous_genotype(bcf_fmt_t *fmt, int gt_idx, int sample_idx, int value) {

    // Generate a pointer to the sample's data in the VCF record
    u_int8_t *ptr = (fmt->p + sample_idx * fmt->size);

    // Generalize the notion of conversion and value writing for various integer sizes
    #define WRITE_HOME_SWITCH_PATH(converter, type) converter(value, (ptr)); converter(value, (ptr + sizeof(type)));

    // Check which integer size is used for the variant, and convert data accordingly
    switch (fmt[gt_idx].type) {
        case BCF_BT_INT8:
            // In the case of an 8-bit (1-byte) integer, Endianness does not matter
            *(ptr) = value;
            *(ptr + sizeof(int8_t)) = value;
            break;
        case BCF_BT_INT16: WRITE_HOME_SWITCH_PATH(i16_to_le, int16_t); break;
        case BCF_BT_INT32: WRITE_HOME_SWITCH_PATH(i32_to_le, int32_t); break;
    }
}

/**
 * The entry point of the program.
 * 
 * @param argc The number of arguments specified in the command-line invocation of the program.
 * @param argv An array of length `argc`, containing pointers to argument strings.
 * @return EXIT_SUCCESS on completion of the program, and EXIT_FAILURE if a runtime error occurs.
 */
int main(int argc, char *argv[]) {

    // Print the proper usage if an incorrect number of arguments is provided
    if (argc != N_ARGS) {
        fprintf(stderr, "Incorrect number of arguments provided!\n");
        fprintf(stderr, "Expected usage:\n\n");
        fprintf(stderr, "refreplacer merged_vcf_path sample_directory_path threads chrom output_vcf_path\n");
        fprintf(stderr, "\tmerged_vcf_path: The path to the merged sample VCF with alternate calls.\n");
        fprintf(stderr, "\tsample_directory_path: The path to the directory containing the sample VCF files.\n");
        fprintf(stderr, "\tthreads: The number of threads.\n");
        fprintf(stderr, "\tchrom: The chromosome or sequence name.\n");
        fprintf(stderr, "\toutput_vcf_path: The path to an output VCF file to create.\n\n");
        exit(EXIT_FAILURE);
    }

    // Define pointers to each argument
    char *merged_vcf_path = argv[ARG_MERGED_VCF_PATH];
    char *sample_directory_path = argv[ARG_SAMPLE_DIRECTORY_PATH];
    char *threads = argv[ARG_THREADS];
    char *chrom = argv[ARG_CHROM];
    char *output_vcf_path = argv[ARG_OUTPUT_VCF_PATH];

    int n_threads = atoi(threads);

    // Read header from BIC's merged file
    htsFile *merged_in = hts_open(merged_vcf_path, "r");
    bcf_hdr_t *out_header = bcf_hdr_read(merged_in);
    hts_close(merged_in);

    // Open output file and write header to it
    htsFile *out = hts_open(output_vcf_path, "wz");
    if (vcf_hdr_write(out, out_header) < 0) {
        fprintf(stderr, "Unable to write to VCF file.\n");
        exit(EXIT_FAILURE);
    }

    // Open a sync reader for multiple VCF files
    bcf_srs_t *sr = bcf_sr_init();
    initialize_sync_reader(sr, merged_vcf_path, sample_directory_path, n_threads, chrom);

    // Create map from sync reader's files to sample header indices in output file
    int sample_indices[sr->nreaders];
    map_readers_to_samples(sr, sample_indices, out_header);

    // Iterate through sync reader until all synced lines are read
    while (bcf_sr_next_line(sr)) {

        if (!bcf_sr_has_line(sr, 0)) {
            // No line present in merged VCF file
            continue;
        }

        // Get output line for merged VCF file
        // Unpack details for all the samples
        bcf1_t *out_line = bcf_sr_get_line(sr, 0);
        bcf_unpack(out_line, BCF_UN_ALL);

        // Get index for GT field in the output line's format field
        bcf_fmt_t *out_fmt = out_line->d.fmt;
        int out_gt_idx = get_gt_index(out_line, out_fmt, out_header);

        if (out_gt_idx < 0) {
            // No GT entry in output line
            continue;
        }

        // Iterate through each sample VCF's reader
        for (unsigned int i = 1; i < sr->nreaders; i++) {

            if (!bcf_sr_has_line(sr, i)) {
                // No line present in sample VCF file
                continue;
            }

            // Get input line from sample VCF file
            // Unpack details for the sample
            bcf1_t *line = bcf_sr_get_line(sr, i);
            bcf_unpack(line, BCF_UN_ALL);

            // Get pointer to sample VCF header
            bcf_hdr_t *sample_header = sr->readers[i].header;
            
            if (strcmp(sample_header->id[BCF_DT_ID][*(line->d.flt)].key, "LowQual") == 0) {
                // Low quality
                continue;
            }

            // Get index for GT field in the input line's format field
            bcf_fmt_t *fmt = line->d.fmt;
            int gt_idx = get_gt_index(line, fmt, sample_header);

            if (gt_idx < 0) {
                // Genotype not present
                continue;
            }

            // Make sure that:
            //  1. The input sample is reference homozygous (0/0) for the current variant and;
            //  2. The merged VCF file has a missing genotype (./.) for the current variant
            if (is_homozygous_for(fmt, gt_idx, 0, REF_ALLELE) && is_homozygous_for(out_fmt, out_gt_idx, sample_indices[i], MISSING_ALLELE)) {

                write_homozygous_genotype(out_fmt, out_gt_idx, sample_indices[i], REF_ALLELE);
            }
        }

        // Write the current line to the output file after inserting all reference genotypes
        if (vcf_write(out, out_header, out_line) < 0) {
            fprintf(stderr, "Unable to write to VCF file.\n");
            exit(EXIT_FAILURE);
        }
    }

    // Destroy thread memory
    if (n_threads > 1) {
        bcf_sr_destroy_threads(sr);
    }

    // Destroy sync reader
    bcf_sr_destroy(sr);

    // Destroy output file
    hts_close(out);

    // Destroy VCF header
    bcf_hdr_destroy(out_header);

    return EXIT_SUCCESS;
}