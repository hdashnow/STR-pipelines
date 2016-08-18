// Bpipe pipeline to run LobSTR allelotyper on BWA bam files
// Set up to run on MCRI Bioinformatics server

/////////////////////////////
// LobSTR

// Paths relative to the analysis directory (or absolute)
LOBSTR_MODELS="/mnt/storage/harrietd/src/lobSTR-bin-Linux-x86_64-4.0.0/share/lobSTR/models"
LOBSTR_BUNDLE="/mnt/storage/harrietd/ref-data/hg19_v3.0.2"
LOBSTR_REF="$LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_ref/lobSTR_"
LOBSTR_STRINFO="$LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_strinfo.tab"
REF='/mnt/storage/shared/genomes/hg19/gatk/bwamem/gatk.ucsc.hg19.fasta'

//load "/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/bpipe.config"

@transform("vcf")
allelotype = {
    transform('.bam') to('.lobstr.vcf', '.lobstr.allelotype.stats') {
        def outname = output.prefix[0..-2]
        exec """
            allelotype \
              --command classify \
              --bam $input.bam \
              --noise_model $LOBSTR_MODELS/illumina_v2.0.3 \
              --out $output1.prefix \
              --strinfo $LOBSTR_STRINFO \
              --index-prefix $LOBSTR_REF
        """, "medium"
    }
}

// -O v for vcf, -O z for vcf.gz
// Also, sort.
@filter('trimmed')
trim_variants = {
        exec "bcftools norm -f $REF -O v $input.vcf | vcf-sort > $output.vcf"
}

run {
    "%.bam" * [
        allelotype + trim_variants
    ]
}
