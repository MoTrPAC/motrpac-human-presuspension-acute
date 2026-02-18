gsutil cp \
    ${snakemake_params[bucket_dir]}/Sample_${snakemake_wildcards[individual]}/analysis/*.final.bam* \
    $TMPDIR/

bcftools mpileup \
    --fasta-ref ${snakemake_params[genome]} \
    --targets-file ${snakemake_input[vcf]} \
    --output-type u \
    --output $TMPDIR/${snakemake_wildcards[individual]}.mpileup \
    $TMPDIR/${snakemake_wildcards[individual]}.final.bam

bcftools call \
    --multiallelic-caller \
    --constrain alleles \
    --targets-file ${snakemake_input[vcf]} \
    --insert-missed \
    --output-type u \
    --output $TMPDIR/${snakemake_wildcards[individual]}.calls.vcf \
    $TMPDIR/${snakemake_wildcards[individual]}.mpileup

bcftools filter \
    --soft-filter LowQual \
    --exclude 'QUAL<20 || DP>100' \
    --output-type z \
    --output ${snakemake_output} \
    $TMPDIR/${snakemake_wildcards[individual]}.calls.vcf
