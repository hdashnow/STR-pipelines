// Paths relative to the analysis directory
PATH_TO_LOBSTR="/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/lobSTR-3.0.3"
LOBSTR_BUNDLE="/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/hg19_v3.0.2"
LOBSTR_REF_PATH="$LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_ref"

load "/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/bpipe.config"

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    branch.lane = 001
    }

@transform("vcf")
allelotype = {
    exec """
        allelotype \
          --command classify \
          --bam $input.bam \
          --noise_model $PATH_TO_LOBSTR/models/illumina_v3.pcrfree \
          --out $output.prefix \
          --strinfo $LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_strinfo.tab \
          --index-prefix $LOBSTR_REF_PATH/lobSTR_
"""
}

run {
    "%.bam" * [
        allelotype
    ]
}
