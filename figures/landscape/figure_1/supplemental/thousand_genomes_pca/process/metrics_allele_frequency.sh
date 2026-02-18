vcftools \
    --gzvcf ${snakemake_input[vcf]} \
    --freq2 \
    --min-alleles 2 \
    --max-alleles 2 \
    --out FILTERED/merged.annotated.missingness_0_02.maf_0_025.hwe_1e-6.snps.allele_frequencies

# Update header to have a column name for reference and alternate frequency
sed -i '1s/{FREQ}/REF_FREQ\tALT_FREQ/g' ${snakemake_output}

# Calculate minor allele frequency
awk -F '\t' '
    BEGIN { OFS="\t"; }
    NR == 1 { print $0, "MAF"; }
    NR > 1 {
        if ($5 < $6) { print $0, $5; }
        else { print $0, $6; } 
    }' ${snakemake_output} > FILTERED/tmp_file
mv FILTERED/tmp_file ${snakemake_output}
