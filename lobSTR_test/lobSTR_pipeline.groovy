// Paths relative to the analysis directory
PATH_TO_LOBSTR="../lobSTR-3.0.3"
LOBSTR_BUNDLE="../hg19_v3.0.2"
LOBSTR_REF_PATH="$LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_ref"

load "../bpipe.config"

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    branch.lane = (input1 =~ /.*L([0-9]*)_*R.*/)[0][1].toInteger()
    }

//@transform("bam")
// Note the replaceAll with this regex will remove the last extension of the filename
lobSTR = {
    from("fastq.gz","fastq.gz") transform(input.prefix.prefix + ".aligned.bam") {
        exec """
            lobSTR \
              --p1 $input1 \
              --p2 $input2 \
              -q --gzip \
              --index-prefix $LOBSTR_REF_PATH/lobSTR_ \
              -o ${output.prefix.replaceAll('\\.[^.]*$','')} \
              --rg-sample ${branch.sample} --rg-lib ${branch.lane}
        """
    }
}

@filter("sort")
sort_bam = {
    from("aligned.bam") {
        exec "samtools sort $input1 $output.prefix"
    }
}

//Not working, and I'm not sure why
merge_bams = {
    doc "Merge BAM files from multiple lanes or samples together. BAM files should have unique sample names and / or read groups"
    produce(branch.sample + "merged.bam") {
        exec """
                MergeSamFiles
                    ${inputs.bam.split().collect { "INPUT="+it }.join(' ')} \
                    USE_THREADING=true  \
                    VALIDATION_STRINGENCY=LENIENT  \
                    AS=true  \
                    OUTPUT=$output
        """
    }
}

index_bam = {
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
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

//run { 
//    "%_R*.fastq.gz" * [ lobSTR + sort_bam + index_bam + index_bam + allelotype 
//}

run {
    ~"(MUW[0-9]*)_.*.fastq.gz" * [
        set_sample_info +
    ~"(L[0-9]*)" * [ lobSTR + sort_bam + index_bam ]
    + merge_bams + index_bam + allelotype
    ]
}
