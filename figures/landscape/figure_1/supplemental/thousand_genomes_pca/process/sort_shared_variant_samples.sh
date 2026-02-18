bcftools view \
    --header \
    ${snakemake_input} | \
    awk '1;/^##INFO=<ID=DP.+$/ {
        print "##INFO=<ID=I16,Number=16,Type=Integer,Description=\"Pileup Information\">";
    }' > \
    $TMPDIR/${snakemake_wildcards[individual]}.unsorted.header.txt

bcftools reheader \
    --header $TMPDIR/${snakemake_wildcards[individual]}.unsorted.header.txt \
    --output $TMPDIR/${snakemake_wildcards[individual]}.unsorted.fixed.vcf.gz \
    ${snakemake_input}

bcftools sort \
    --max-mem 15.9G \
    --temp-dir $TMPDIR/ \
    --output-type z \
    --output ${snakemake_output} \
    $TMPDIR/${snakemake_wildcards[individual]}.unsorted.fixed.vcf.gz
