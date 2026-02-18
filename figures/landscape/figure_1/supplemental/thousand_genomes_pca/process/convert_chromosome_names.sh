bcftools view ${snakemake_input} | \
    awk '{ if ($0 ~ /^MT/) { print "chrM"; } else if ($0 !~ /^#/) { print "chr"$0; } else { print $0; } }' | \
    bcftools view \
    --threads ${snakemake[threads]} \
    --output-type z \
    --output ${snakemake_output} \
    -