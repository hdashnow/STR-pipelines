// Bpipe pipeline to simulate paired end reads from a fasta file, with no stutter
// Simulates reads from STR alleles at a provided set of loci

// Set up to run on the redwood cluster

INSTALLDIR="~/storage/git/STR-pipelines"
ART="~/tools/art_bin_MountRainier"
REF="/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/ref-data/GATK_Bundle/human_g1k_v37_decoy.fasta"
TOOLS="$INSTALLDIR/simulate_reads" //custom R/python scripts
STUTTER="$INSTALLDIR/simulate_reads/no_stutter_model.csv"
GATK="~/tools/gatk-4.1.2.0/gatk"

DECOY_REF="$INSTALLDIR/simulate_reads/reference-data/hg19.STRdecoys.sorted.fasta"
EXOME_TARGET="/group/bioi1/harrietd/ref-data/hg19_RefSeq_coding.sorted.bed"

// Adjust simulation parameters
PLATFORM="illumina"
total_coverage = 30

def get_fname(path) {
    def x = path.split("/")[-1]
    return(x)
}

@preserve("*.bed")
mutate_locus = {
    doc """Generate a random heterozygous coding STR loci in the potentially
        pathogenic range."""

    output.dir = "sim_bed"
    branch.simID = branch.name

    produce(branch.simID + ".bed") {
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
            bedtools sort -i $input.bed > $output.bed
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
                python $TOOLS/generate_stutter_vcfs.py $REF $input.bed --output $output.prefix.prefix --stutter $STUTTER --flank 100 > $output.txt
        """
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

index_vcf = {
    doc "produces .vcf.idx"
    output.dir = "vcf_bed"

produce(input.vcf + ".idx") {
        exec """
            $GATK IndexFeatureFile
                -F $input.vcf
        """
}
    forward input
}

@Transform("fasta")
mutate_ref = {
    doc "Generate a version of the reference genome (or subset) with mutations given by the input VCF"
    output.dir = "fasta"

    exec """
        $GATK FastaAlternateReferenceMaker
            -R $REF
            -O $output.fasta
            -L $input.bed
            -V $input.vcf
    """
}


/////////////////////////////
// Generate reads

generate_reads = {
    doc "Sample reads from the altered reference sequence segment"
    output.dir = "fastq"

    produce( get_fname(input.fasta.prefix) + "_L001_R1.fq", get_fname(input.fasta.prefix) + "_L001_R2.fq") {

        // Set target coverage for this stutter allele
        def coverage = 0.5 * total_coverage
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
    output.dir = "fastq"

    produce(input1.prefix + ".fastq.gz") {
        exec "cat $inputs.fq | gzip -c > $output.gz"
    }
}


/////////////////////////////
// Align reads


@preserve("*.bam")
align_bwa = {
    doc "Concatenate with background reads then align with bwa mem algorithm."

    def fname = get_fname(input1)
    def lane = "001"
    def sample = branch.name
    from("fastq.gz", "fastq.gz") produce(fname.prefix.prefix + ".bam") {
        exec """
            bwa mem -M
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample}\\tSM:${sample}"
            $REF
            $input1.gz
            $input2.gz |
            samtools view -bSuh - | samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwamem"
    }
}

index_bam = {
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

// -O v for vcf, -O z for vcf.gz
// Also, sort.
@filter("trimmed")
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
            merge_bed + index_vcf + mutate_ref + generate_reads
        ] +

        "%0.5_*.stutter.merged_L001_R%.fq" * [
            combine_gzip
        ] +

        "%_R*.fastq.gz" * [
            align_bwa + index_bam
        ] +

        "%.truth.vcf" * [ trim_variants ] 
}
