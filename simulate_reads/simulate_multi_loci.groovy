// Bpipe pipeline to simulate paired end reads from a fasta file, with no stutter
// Simulates reads from STR alleles at a provided set of loci

// Set up to run on the redwood cluster

INSTALLDIR='~/storage/git/STR-pipelines'
ART='~/tools/art_bin_MountRainier'
REF='/group/bioi1/shared/genomes/hg19/gatk/bwamem/gatk.ucsc.hg19.fasta'
CHR_ORDER='/group/bioi1/shared/genomes/hg19/gatk/gatk.ucsc.hg19.chr_order.txt'
TOOLS='$INSTALLDIR/simulate_reads' //custom R/python scripts
STUTTER='$INSTALLDIR/simulate_reads/no_stutter_model.csv'

DECOY_REF="$INSTALLDIR/simulate_reads/reference-data/hg19.STRdecoys.sorted.fasta"
EXOME_TARGET="/group/bioi1/harrietd/ref-data/hg19_RefSeq_coding.sorted.bed"

// Adjust simulation parameters
PLATFORM='illumina'
total_coverage = 30
//LOCUS='chr2:233712201-233712246'
LOCUS='chr13:70713515-70713561' //ATXN80OS CAG repeat
// Adjust number of variants to simulate here
simID = (1..100)

def get_fname(path) {
    x = path.split("/")[-1]
    return(x)
}

param_map = [:] //Define outside function so one big map for all branches

// Parses tab-delmited table of filenames and their corresponding parameters
// Returns these as a map that can be used to get the correct parameter as a branch variable
def parse_parameters(all_parameters) {

    def lines = all_parameters.readLines()
    lines.each {
        def row_list = it.split()
        param_map[get_fname(row_list[0])] = ['probability':row_list[1], 'bedfile':row_list[2]]
    }
}

@preserve("*.bed")
mutate_locus = {
    doc """Generate a random heterozygous coding STR loci in the potentially
        pathogenic range."""

    output.dir = "sim_bed"
    branch.simID = branch.name

    produce(branch.simID + '.bed') {
        exec """
            $INSTALLDIR/simulate_reads/STR_simulation_script.R
                -L $LOCUS
                $INSTALLDIR/simulate_reads/reference-data/hg19.simpleRepeat.txt.gz
                $INSTALLDIR/simulate_reads/reference-data/str-stats
                -O $output.bed
                -m 500
        """

    }
}

@filter("sorted")
sort_bed = {
    doc "sort bed file"
    output.dir = "sim_bed"
    branch.source_bed = input.bed

    preserve("*.bed") {
        exec """
            bedtools sort -i $input.bed -faidx $CHR_ORDER > $output.bed
        """
    }
}



/////////////////////////////
// Produce mutated fasta file

generate_vcf = {
    doc "Generate a VCF of STR mutations and stutter, along with their probabilities"
    output.dir = "vcf_bed"

    def bedname = get_fname(input.bed)

    preserve("*.truth.vcf") {
        produce(bedname.prefix + ".truth.vcf", "*.stutter.vcf", bedname.prefix + ".txt", "*.stutter.bed") {

            exec """
                /group/bioi1/harrietd/src/miniconda3/envs/STR/bin/python $TOOLS/generate_stutter_vcfs.py $REF $input.bed --output $output.prefix.prefix --stutter $STUTTER > $output.txt
        """

        // Add settings for this branch param_map
        new File("$output.dir").listFiles()
        def File all_params = new File( output.txt )
        parse_parameters(all_params)

        }
    }
}

@filter("merged")
merge_bed = {
    doc "merge bed file"
    output.dir = "vcf_bed"

    exec """
        bedtools merge -i $input.bed > $output.bed
    """

//    forward input
}

@Transform("fasta")
mutate_ref = {
    doc "Generate a version of the reference genome (or subset) with mutations given by the input VCF"
    output.dir = "fasta"

    exec """
        java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar
            -T FastaAlternateReferenceMaker
            -R $REF
            -o $output.fasta
            -L $input.bed
            -V $input.vcf
    """
}


/////////////////////////////
// Generate reads

generate_reads = {
    doc "Sample reads from the altered reference sequence segment"
    output.dir = "fastq"

    produce( get_fname(input.fasta.prefix) + '_L001_R1.fq', get_fname(input.fasta.prefix) + '_L001_R2.fq') {

        // Set target coverage for this stutter allele
//println("$input.vcf $param_map")
        def coverage = param_map[get_fname("$input.vcf")]["probability"].toDouble() * total_coverage
        def outname = output.prefix[0..-2]
        exec """
            $ART/art_illumina -i $input.fasta -p -na
                -l 150 -ss HS25 -f $coverage
                -m 500 -s 50
                -o $outname
        """
    }
}

combine_gzip = {
    def ID = get_fname(input1).split("\\.")[0]
    output.dir = "fastq"

    from('*.fq') produce(ID + input1.prefix[-8..-1] + '.fastq.gz') {
        preserve("*.gz") {
            exec "cat $inputs.fq | gzip -c > $output.gz"
        }
    }
}


/////////////////////////////
// Align reads

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    branch.lane = 001
    }

threads=8

@preserve("*.bam")
align_bwa = {
    doc "Concatenate with background reads then align with bwa mem algorithm."

    def fastaname = get_fname(REF)
    from('fastq.gz', 'fastq.gz') produce(branch.name + '.bam') {
        exec """
            bwa mem -M -t $threads
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample}\\tSM:${sample}"
            $DECOY_REF
            $input1.gz
            $input2.gz |
            samtools view -bSuh - | samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwamem"
    }
}

@preserve("*.bai")
index_bam = {
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

// -O v for vcf, -O z for vcf.gz
// Also, sort.
@filter('trimmed')
trim_variants = {
    preserve("*.vcf") {
        exec "bcftools norm -f $REF -O v $input.vcf | vcf-sort > $output.vcf"
    }
}

clean_intermediates = { cleanup "*.fq", "*.fastq.gz", "*.stutter.vcf", "*.stutter.bed" "*.stutter.merged.fasta"}

/////////////////////////////
// Run pipeline

run {
        sort_bed +
        generate_vcf +

        "%.stutter.*" * [
            merge_bed + mutate_ref + generate_reads
        ] +

        "%.sorted.*.stutter.merged_L001_R%.fq" * [
            combine_gzip
        ] +

        '%_R*.fastq.gz' * [
            set_sample_info +
            align_bwa + index_bam
        ] +

        "%.truth.vcf" * [ trim_variants ] //+
        //clean_intermediates
}
