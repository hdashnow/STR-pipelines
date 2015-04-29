PATH_TO_LOBSTR="/vlsci/VR0320/hdashnow/test_msat_genotyping/lobSTR_test/lobSTR-3.0.3"
LOBSTR_BUNDLE="/vlsci/VR0320/hdashnow/test_msat_genotyping/lobSTR_test/hg19_v3.0.2"
LOBSTR_REF_PATH="$LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_ref"

load "/vlsci/VR0320/hdashnow/test_msat_genotyping/lobSTR_test/bpipe.config"

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    branch.lane = (input1 =~ /.*L([0-9]*)_*R.*/)[0][1].toInteger()
    }

//@transform("bam")
// Note the replaceAll with this regex will remove the last extension of the filename
map = {
    from("fastq.gz","fastq.gz") transform(input.prefix.prefix + ".aligned.bam") {
        exec """
            echo input1 $input1 \n
            echo input2 $input2 \n
            echo output ${output.prefix.replaceAll('\\.[^.]*$','')} \n
            echo $inputs > ${output.prefix.replaceAll('\\.[^.]*$','')}.aligned.bam \n
            echo rg-sample ${branch.sample} --rg-lib ${branch.lane} \n
        """
    }
}

//Not working, and I'm not sure why
@filter("merge")
merge_bams = {
    produce(branch.sample + ".bam") {
        doc "Merge BAM files from multiple lanes or samples together. BAM files should have unique sample names and / or read groups"
        exec """
            cat $inputs.bam > ${branch.sample}.bam
        """
    }
}

@transform("vcf")
allelotype = {
    exec """
    cat $inputs.bam > $output.vcf     
"""
}

run { 
    //"%_R*.fastq.gz" * [ 
    ~"(MUW[0-9]*)_.*.fastq.gz" * [
        set_sample_info +
    ~"(L[0-9]*)" * [ map ]
    + merge_bams + allelotype
    ]
}
