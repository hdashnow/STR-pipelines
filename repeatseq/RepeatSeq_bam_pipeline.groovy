// Paths relative to the analysis directory
load "/vlsci/VR0320/hdashnow/git/STR-pipelines/repeatseq/bpipe.config"
// RepeatSeq installation
REPEATSEQ="/vlsci/VR0320/hdashnow/git/repeatseq"

// Reference data
REF="/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta" // Genome
REGIONS="$REPEATSEQ/regions/hg19.2014.regions"
// STRs in the format required by RepeatSeq

PLATFORM='illumina'

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    branch.lane = (input1 =~ /.*L([0-9]*)_*R.*/)[0][1].toInteger()
}

align_bwa = {
    doc "Align with bwa mem algorithm."
    from('fastq.gz', 'fastq.gz') transform('bam') {
        exec """
            bwa mem -M -t $threads
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample}\\tSM:${sample}"  
            $REF $input1.gz $input2.gz | 
            samtools view -bSuh - | samtools sort - $output.prefix
        """, "bwamem"
    }
}

index_bam = {
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

genotype = {

    doc "Genotype STRs using RepeatSeq"

    from("bam") transform("bam.vcf") {
        exec """
            $REPEATSEQ/repeatseq $input.bam $REF $REGIONS
        """
    }
}

run {
//    '%_R*.fastq.gz' * [
//        set_sample_info +
//        align_bwa + index_bam
//    ] +

    "%.bam" * [
        genotype
    ]
}
