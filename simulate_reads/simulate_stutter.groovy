// Bpipe pipeline to simulate paired end reads from a fasta file, with microsatellite stutter

ART='/vlsci/VR0320/hdashnow/art_bin_ChocolateCherryCake'
REF='/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta' 
GATK="/usr/local/gatk/3.4-46/GenomeAnalysisTK.jar"

def get_fname(path) {
    x = path.split("/")[-1]
    return(x)
}

/////////////////////////////
// Produce mutated fasta file

//XXX TO DO
generate_vcf = {
    doc "Generate a VCF of STR mutations and stutter"
}

@Transform("fasta")
mutate_ref = {
    doc "Generate a version of the reference genome (or subset) with mutations given by the input VCF"
    exec """
        java -Xmx4g -jar $GATK
            -T FastaAlternateReferenceMaker
            -R $REF
            -o $output.fasta
            -L chr6:16,243,709-16,793,549
            -V $input.vcf
    """
}


/////////////////////////////
// Generate reads

generate_reads = {
    def samplename = branch.name
    branch.coverage = branch.name.toInteger()
    def fastaname = get_fname(REF)
    produce(fastaname + '_cov' + samplename + '_L001_R1.fq', fastaname + '_cov' + samplename + '_L001_R2.fq') {
        def outname = output.prefix[0..-2]
        exec """
            $ART/art_illumina -i $input.fasta -p -na 
                -l 150 -ss HS25 -f $coverage
                -m 500 -s 50 
                -o $outname
        ""","medium"
    }
}

gzip = {
    transform('.fq') to('.fastq.gz') {
        exec "gzip -c $input.fq > $output.gz","medium" 
    }
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip') {
        exec "fastqc -o ${output.dir} $inputs.gz","medium"
    }
}

/////////////////////////////
// Align reads

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

// Paths relative to the analysis directory
load "/vlsci/VR0320/hdashnow/git/STR-pipelines/repeatseq/bpipe.config"

PLATFORM='illumina'

threads=8

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

/////////////////////////////
// Run pipeline

coverages = [10]

run {

    coverages * [
        mutate_ref + generate_reads
    ] +

    "%.fq" * [ 
        gzip 
    ] +

    '%_R*.fastq.gz' * [
        set_sample_info +
        align_bwa + index_bam +
        RealignerTargetCreator + IndelRealigner
    ]
}
