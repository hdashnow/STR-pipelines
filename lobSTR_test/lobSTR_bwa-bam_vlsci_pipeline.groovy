// Bpipe pipeline to run LobSTR allelotyper on BWA bam files
// Set up to run on Barcoo, VLSCI

/////////////////////////////
// LobSTR

// Paths relative to the analysis directory (or absolute)
LOBSTR_CODE="/vlsci/VR0320/hdashnow/software/lobSTR-4.0.0"
LOBSTR_BUNDLE="/vlsci/VR0320/storage/hdashnow/ref-data/hg19_v3.0.2"
LOBSTR_REF="$LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_ref/lobSTR_"
LOBSTR_STRINFO="$LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_strinfo.tab"

//load "/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/bpipe.config"

@transform("vcf")
allelotype = {
    transform('.bam') to('.lobstr.vcf', '.lobstr.allelotype.stats') {
        def outname = output.prefix[0..-2]
        exec """
            allelotype \
              --command classify \
              --bam $input.bam \
              --noise_model $LOBSTR_CODE/models/illumina_v2.0.3 \
              --out $output1.prefix \
              --strinfo $LOBSTR_STRINFO \
              --index-prefix $LOBSTR_REF
        """, "medium"
    }
}


run {
    "%.bam" * [
        allelotype
    ]
}

