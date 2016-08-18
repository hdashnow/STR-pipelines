// Bpipe pipeline to simulate paired end reads from a fasta file

ART='/vlsci/VR0320/hdashnow/art_bin_ChocolateCherryCake'
REF='/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta' 
GATK='/usr/local/gatk/3.4-46/GenomeAnalysisTK.jar'

generate_reads = {
    def samplename = branch.name
    produce(samplename + '_L001_R1.fq', samplename + '_L001_R2.fq') {
        def outname = output.prefix[0..-2]
        exec """
            $ART/art_illumina -i $input.fasta -p -na 
                -l 150 -ss HS25 -f 30 
                -m 500 -s 50 
                -o $outname
        ""","quick"
    }
}

gzip = {
    transform('.fq') to('.fastq.gz') {
        exec "gzip -c $input.fq > $output.gz","quick" 
    }
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip') {
        exec "fastqc -o ${output.dir} $inputs.gz","quick"
    }
}

PLATFORM='illumina'

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    branch.lane = 001 //(input1 =~ /.*L([0-9]*)_*R.*/)[0][1].toInteger()
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

RealignerTargetCreator = {
    transform("bam") to ("intervals") {
        exec """
            java -Xmx2g -jar $GATK 
                -T RealignerTargetCreator
                -R $REF
                -I $input.bam
                -o $output.intervals
        """
    }
    forward input
}

@filter("realigned")
IndelRealigner = {
    exec """
        java -Xmx4g -jar $GATK
            -T IndelRealigner
            -R $REF
            -I $input.bam
            -targetIntervals $input.intervals
            -o $output.bam
            --consensusDeterminationModel USE_SW
            --maxPositionalMoveAllowed 500
    """
}


run {
    "%.fasta" * [
        generate_reads
    ] +
    "%.fq" * [ 
        gzip //+ fastqc
    ] +
    '%_R*.fastq.gz' * [
        set_sample_info + 
        align_bwa + index_bam +
        RealignerTargetCreator + IndelRealigner
    ]
}
