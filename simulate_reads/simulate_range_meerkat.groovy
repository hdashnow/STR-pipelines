// Bpipe pipeline to simulate paired end reads from a fasta file, with no microsatellite stutter
// Simulates long STRs at a selected set of loci, and
// simulates normal STRs over the rest of the genome (background)
// Variant STRs simulated over a range of sizes, repeat units (and coverage?)
// Set up to run on Meerkat cluster at MCRI

ART='/group/bioi1/harrietd/src/art_bin_MountRainier'
REF='/group/bioi1/shared/genomes/hg19/gatk/bwamem/gatk.ucsc.hg19.fasta'
CHR_ORDER='/group/bioi1/shared/genomes/hg19/gatk/gatk.ucsc.hg19.chr_order.txt'
TOOLS='/group/bioi1/harrietd/git/STR-pipelines/simulate_reads'
STUTTER='/group/bioi1/harrietd/git/STR-pipelines/simulate_reads/no_stutter_model.csv'

DECOY_REF="/group/bioi1/harrietd/git/STR-pipelines/simulate_reads/reference-data/hg19.STRdecoys.fasta"
EXOME_TARGET="/group/bioi1/harrietd/ref-data/hg19_RefSeq_coding.bed"

PLATFORM='illumina'
total_coverage = 50

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
mutate_all = {
    doc """Generate a random heterozygous coding STR loci in the potentially
        pathogenic range."""

    output.dir = "sim_bed"
    branch.simID = branch.name

    produce(branch.simID + '.mutant.bed', branch.simID + '.background.bed') {
        exec """
            /group/bioi1/harrietd/git/STR-pipelines/simulate_reads/STR_simulation_script.R
                -L random
                /group/bioi1/harrietd/git/STR-pipelines/simulate_reads/reference-data/hg19.simpleRepeat.txt.gz
                /group/bioi1/harrietd/git/STR-pipelines/simulate_reads/reference-data/str-stats
                -O $output1.bed
                --background $output2.bed
        """

    }
}

cat_bed = {
    doc "combine multiple bed files and sort them"
    output.dir = "sim_bed"

    produce(branch.simID + '.all.bed') {
        exec """
            cat $inputs.bed | bedtools sort -faidx $CHR_ORDER -i stdin  > $output.bed
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
        cut -f 1,2,3,4 $EXOME_TARGET | cat - $input.bed | bedtools sort -i stdin | bedtools merge -i stdin > $output.bed
    """
    
    forward input
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

@preserve("*.txt")
STR_coverage = {
    transform("bam") to ("coverage.txt") {
        exec """
            bedtools coverage -counts
            -a /group/bioi1/harrietd/git/STR-pipelines/simulate_reads/reference-data/STRdecoys.bed
            -b $input.bam > $output.txt
        """
    }
}

STR_locus_counts = {
    transform("bam") to ("locus_counts.txt") {
        exec """
            /group/bioi1/harrietd/src/miniconda3/envs/STR/bin/python /group/bioi1/harrietd/git/micro-genotyper-long/python_code/identify_locus.py
            --bam $input.bam
            --bed /group/bioi1/harrietd/git/STR-pipelines/simulate_reads/reference-data/hg19.simpleRepeat_period1-6.bed
            --output $output.txt
        """
    }
}

@transform('coverage')
coverage = {
    exec """
        bedtools coverage
            -sorted -d 
            -g ${REF}.genome
            -a $EXOME_TARGET 
            -b $input.bam 
            > $output.coverage
     ""","bedtools"
}

median_cov = {
    exec """
        cut -f 5 $input.coverage  | sort -n | awk -f /group/bioi1/harrietd/git/micro-genotyper-long/repeat_genotyper_bpipes/median.awk > $output
     """
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

// Adjust number of variants to simulate here
simID = (1..5)

run {
// Generate bed file of loci to simulate. One locus pathogenic per file (rest normal range).
// Do this multiple times until a few scenarios have been covered
// Also one bed file of all the remaining loci (not in the above file) with normal ranges (background).

    simID * [

        mutate_all  +
        cat_bed +
        generate_vcf 

    ] +

        "%.stutter.*" * [
            merge_bed + mutate_ref + generate_reads
        ] +

        "%.all.*.stutter_L001_R%.fq" * [
            combine_gzip
        ] +

        '%_R*.fastq.gz' * [
            set_sample_info +
            align_bwa + index_bam +
            STR_locus_counts +
            coverage +
            median_cov +
            STR_coverage
        ] +

        "%.truth.vcf" * [ trim_variants ] //+
        //clean_intermediates
}
