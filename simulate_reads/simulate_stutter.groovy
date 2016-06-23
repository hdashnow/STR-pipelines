// Bpipe pipeline to simulate paired end reads from a fasta file, with microsatellite stutter

ART='/vlsci/VR0320/hdashnow/art_bin_ChocolateCherryCake'
REF='/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta'
GATK="/usr/local/gatk/3.4-46/GenomeAnalysisTK.jar"
TOOLS='/vlsci/VR0320/hdashnow/git/STR-pipelines/simulate_reads'
STUTTER='/vlsci/VR0320/hdashnow/git/STR-pipelines/simulate_reads/stutter_model.csv'
total_coverage = 100

def get_fname(path) {
    x = path.split("/")[-1]
    return(x)
}

// Parses tab-delmited table of filenames and their corresponding parameters
// Returns these as a map that can be used to get the correct parameter as a branch variable
def parse_parameters(all_parameters) {

    def param_map = [:]
    lines = all_parameters.readLines()

    lines.each {
        row_list = it.split()
        param_map[row_list[0]] = ['probability':row_list[1], 'bedfile':row_list[2], 'delta':row_list[3]]
    }
    return(param_map)
}

/////////////////////////////
// Produce mutated fasta file

generate_vcf = {
    doc "Generate a VCF of STR mutations and stutter, along with their probabilities"

    def bedname = get_fname(input.bed)

    preserve("*.truth.vcf") {
        produce(bedname.prefix + ".truth.vcf", "*.vcf", bedname.prefix + ".txt") {
            exec """
                python $TOOLS/generate_stutter_vcfs.py $REF $input.bed --output $output.prefix.prefix --stutter $STUTTER > $output.txt
        """
        File all_params = new File( output.txt )
        branch.param_map = parse_parameters(all_params)
        }
    }
}

@Transform("fasta")
mutate_ref = {
    doc "Generate a version of the reference genome (or subset) with mutations given by the input VCF"

        // Set target coverage for this stutter allele
        branch.coverage = branch.param_map["$input.vcf"]["probability"].toDouble() * total_coverage
        branch.bedfile = branch.param_map["$input.vcf"]["bedfile"]
        branch.delta = branch.param_map["$input.vcf"]["delta"]
    exec """
        java -Xmx4g -jar $GATK
            -T FastaAlternateReferenceMaker
            -R $REF
            -o $output.fasta
            -L $branch.bedfile
            -V $input.vcf
    ""","small"
}


/////////////////////////////
// Generate reads

generate_reads = {
    def readname = 'stutter_' + branch.delta + '_'
    def fastaname = get_fname(REF)
    produce(fastaname + '_' + input.fasta.prefix + '_L001_R1.fq', fastaname + '_' + input.fasta.prefix + '_L001_R2.fq') {
        def outname = output.prefix[0..-2]
        exec """
            $ART/art_illumina -i $input.fasta -p -na
                -l 150 -ss HS25 -f $branch.coverage
                -m 500 -s 50
                --id $readname 
                -o $outname
        ""","small"
    }
}

combine_gzip = {
    from('*.fq') produce(input.fq.prefix + '.fastq.gz') {
        exec "cat $inputs.fq | gzip -c > $output.gz","small"
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

@preserve("*.bam")
align_bwa = {
    doc "Align with bwa mem algorithm."
    from('fastq.gz', 'fastq.gz') transform('bam') {
        exec """
            bwa mem -M -t $threads
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample}\\tSM:${sample}"
            $REF $input1.gz $input2.gz |
            samtools view -bSuh - | samtools sort -o $output.bam -
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


run {

    generate_vcf +

    "%.stutter_%.vcf" * [
        mutate_ref + generate_reads
    ] +

    "*.stutter_*_R%.fq" * [
        combine_gzip
    ] +

    '%_R*.fastq.gz' * [
        set_sample_info +
        align_bwa + index_bam //+
//        RealignerTargetCreator + IndelRealigner
    ]
}
