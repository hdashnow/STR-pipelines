// Bpipe pipeline to simulate paired end reads from a fasta file, with microsatellite stutter
// Simulates long STRs at a selected set of loci, and
// simulates normal STRs over the rest of the genome (background)

ART='/mnt/storage/harrietd/src/art_bin_MountRainier'
REF='/mnt/storage/shared/genomes/hg19/gatk/bwamem/gatk.ucsc.hg19.fasta'
CHR_ORDER='/mnt/storage/shared/genomes/hg19/gatk/gatk.ucsc.hg19.chr_order.txt'
GATK='/mnt/storage/harrietd/src/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar'
TOOLS='/mnt/storage/harrietd/git/STR-pipelines/simulate_reads'
STUTTER='/mnt/storage/harrietd/git/STR-pipelines/simulate_reads/no_stutter_model.csv'
total_coverage = 100

def get_fname(path) {
    x = path.split("/")[-1]
    return(x)
}


param_map = [:] //Define outside function so one big map for all branches

// Parses tab-delmited table of filenames and their corresponding parameters
// Returns these as a map that can be used to get the correct parameter as a branch variable
def parse_parameters(all_parameters) {

    lines = all_parameters.readLines()

    lines.each {
        row_list = it.split()
        param_map[row_list[0]] = ['probability':row_list[1], 'bedfile':row_list[2]]
    }
    //return(param_map)
}

@preserve("*.bed")
mutate_background = {
    doc """ Generate a set of coding STR loci in the normal range.
        Exclude one locus that will be simulated as variable."""

    output.dir = "sim_bed"

    produce('background.bed') {
        exec """
            /mnt/storage/harrietd/git/STR-pipelines/simulate_reads/STR_simulation_script.R
                -L chr2:233712201-233712246
                /mnt/storage/harrietd/git/STR-pipelines/simulate_reads/reference-data/hg19.simpleRepeat.txt.gz
                /mnt/storage/harrietd/git/STR-pipelines/simulate_reads/reference-data/str-stats
                --background $output.bed > /dev/null
        """

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
            /mnt/storage/harrietd/git/STR-pipelines/simulate_reads/STR_simulation_script.R
                -L chr2:233712201-233712246
                /mnt/storage/harrietd/git/STR-pipelines/simulate_reads/reference-data/hg19.simpleRepeat.txt.gz
                /mnt/storage/harrietd/git/STR-pipelines/simulate_reads/reference-data/str-stats
                -O $output.bed
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
                python $TOOLS/generate_stutter_vcfs.py $REF $input.bed --output $output.prefix.prefix --stutter $STUTTER > $output.txt
        """
        File all_params = new File( output.txt )
        //branch.param_map = parse_parameters(all_params)
        new File("$output.dir").listFiles()
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
}

@Transform("fasta")
mutate_ref = {
    doc "Generate a version of the reference genome (or subset) with mutations given by the input VCF"
    output.dir = "fasta"

        // Set target coverage for this stutter allele
        //println("$name $input.vcf $param_map")
        //branch.coverage = branch.param_map["$input.vcf"]["probability"].toDouble() * total_coverage
        branch.coverage = param_map["$input.vcf"]["probability"].toDouble() * total_coverage
        //branch.bedfile = branch.param_map["$input.vcf"]["bedfile"]
    exec """
        java -Xmx4g -jar $GATK
            -T FastaAlternateReferenceMaker
            -R $REF
            -o $output.fasta
            -L $input.bed
            -V $input.vcf
    ""","quick"
}


/////////////////////////////
// Generate reads

generate_reads = {
    doc "Sample reads from the altered reference sequence segment"
    output.dir = "fastq"

    produce( get_fname(input.fasta.prefix) + '_L001_R1.fq', get_fname(input.fasta.prefix) + '_L001_R2.fq') {
        def outname = output.prefix[0..-2]
        exec """
            $ART/art_illumina -i $input.fasta -p -na
                -l 150 -ss HS25 -f $branch.coverage
                -m 500 -s 50
                -o $outname
        ""","quick"
    }
}

combine_gzip = {
    from('*.fq') produce(input.fq.prefix + '.fastq.gz') {
        preserve("*.gz") {
            exec "cat $inputs.fq | gzip -c > $output.gz","small"
        }
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

@preserve("*.bai")
index_bam = {
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam", "medium"
    }
    forward input
}


DECOY_REF="/mnt/storage/harrietd/git/STR-pipelines/simulate_reads/reference-data/ucsc.hg19.STRdecoys.fasta"
BACKGROUND_R1="/mnt/storage/harrietd/git/STR-pipelines/simulate_reads/simulate_exome/background_coding_path_STR/fastq/background.sorted.0.5_1.stutter.merged_L001_R1.fastq.gz"
BACKGROUND_R2="/mnt/storage/harrietd/git/STR-pipelines/simulate_reads/simulate_exome/background_coding_path_STR/fastq/background.sorted.0.5_1.stutter.merged_L001_R2.fastq.gz"
PLATFORM='illumina'
threads=8

@preserve("*.bam")
align_bwa = {
    doc "Concatenate with background reads then align with bwa mem algorithm."

    def fastaname = get_fname(REF)
    from('fastq.gz', 'fastq.gz') produce(fastaname.prefix + '.' + get_fname(branch.source_bed).prefix + '.bam') {
        exec """
            bwa mem -M -t $threads
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample}\\tSM:${sample}"
            $DECOY_REF
            <(gzip -dc $BACKGROUND_R1 $input1.gz)
            <(gzip -dc $BACKGROUND_R2 $input2.gz) |
            samtools view -bSuh - | samtools sort -o $output.bam -
        """, "bwamem"
    }
}

@preserve("*.txt")
STR_coverage = {
    transform("bam") to ("coverage.txt") {
        exec """
            bedtools coverage -counts
            -a /mnt/storage/harrietd/git/STR-pipelines/simulate_reads/reference-data/STRdecoys.bed
            -b $input.bam > $output.txt
        """
    }
}



// -O v for vcf, -O z for vcf.gz
// Also, sort.
@filter('trimmed')
trim_variants = {
    preserve("*.vcf") {
        exec "bcftools norm -f $REF -O v $input.vcf | vcf-sort > $output.vcf"
    }
}

/////////////////////////////
// Run pipeline

// Adjust number of variants to simulate here
simID = (1..100)

run {
// Generate bed file of loci to simulate. One locus pathogenic per file (rest normal range).
// Do this multiple times until a few scenarios have been covered
// Also one bed file of all the remaining loci (not in the above file) with normal ranges (background).


    simID * [

        mutate_locus  +

        sort_bed +
        generate_vcf +

        "%.stutter.*" * [
            merge_bed + mutate_ref + generate_reads
        ] +

        "*.stutter.*_R%.fq" * [
            combine_gzip
        ] +

        '%_R*.fastq.gz' * [
            set_sample_info +
            align_bwa + index_bam +
            STR_coverage
        ] +

        "%.truth.vcf" * [ trim_variants ]

    ]
}
