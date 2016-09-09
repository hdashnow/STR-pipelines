// Bpipe pipeline to run LobSTR allelotyper on BWA bam files
// Set up to run on the Broad cluster

REF='/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta' 


def get_fname(path) {
    x = path.split("/")[-1]
    return(x)
}
/////////////////////////////
// Steps in common
set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    branch.lane = 001
}

index_bam = {
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam", "medium"
    }
    forward input
}


/////////////////////////////
// LobSTR

// Paths relative to the analysis directory (or absolute)
LOBSTR_CODE="/home/unix/hdashnow/software/lobSTR-bin-Linux-x86_64-4.0.0/share/lobSTR"
LOBSTR_BUNDLE="/home/unix/hdashnow/storage/ref-data/GRCh37_v3.0.2"
LOBSTR_REF="$LOBSTR_BUNDLE/lobstr_v3.0.2_GRCh37_ref/lobSTR_"
LOBSTR_STRINFO="$LOBSTR_BUNDLE/lobstr_v3.0.2_GRCh37_strinfo.tab"

//load "/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/bpipe.config"


@filter("sort")
sort_bam = {
    from("aligned.bam") {
        exec "samtools sort $input1 $output.prefix", "medium"
    }
}

@transform("vcf")
allelotype_single = {

    doc "Run LobSTR allelotyper STR calling on a single bam file at a time"

    from('.bam') transform('.lobstr.vcf', '.lobstr.allelotype.stats') {
        def outname = output.prefix[0..-2]
        exec """
            allelotype \
              --command classify \
              --bam $input.bam \
              --noise_model $LOBSTR_CODE/models/illumina_v3.pcrfree \
              --out $output1.prefix \
              --strinfo $LOBSTR_STRINFO \
              --index-prefix $LOBSTR_REF
        """, "medium"
    }
}

@transform("vcf")
allelotype_multi = {

    doc "Run LobSTR allelotyper STR calling on multiple bam files at the same time"

    from('*.bam') transform('.lobstr.vcf', '.lobstr.allelotype.stats') {
        def outname = output.prefix[0..-2]
        def infiles = inputs.join(',')
        exec """
            allelotype \
              --command classify \
              --bam $infiles \
              --noise_model $LOBSTR_CODE/models/illumina_v3.pcrfree \
              --out $output1.prefix \
              --strinfo $LOBSTR_STRINFO \
              --index-prefix $LOBSTR_REF
        """, "medium"
    }
}

run {
    "%.bam" * [
        allelotype_single
    ]
}

//run {
//    [
//        allelotype_multi
//    ]
//}



