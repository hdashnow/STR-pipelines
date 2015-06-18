// Paths relative to the analysis directory
load "../bpipe.config"
// RepeatSeq installation
REPEATSEQ="/vlsci/VR0320/hdashnow/git/repeatseq"

// Reference data
REF="/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta" // Genome
REGIONS="$REPEATSEQ/regions/hg19.2014.regions"
// STRs in the format required by RepeatSeq

@transform("vcf")
genotype = {

    doc "Genotype STRs using RepeatSeq"

    produce("$output.vcf") {
        exec """
            $REPEATSEQ/repeatseq $input.bam $REF $REGIONS
        """
    }
}

run {
    "%.bam" * [
        genotype
    ]
}
