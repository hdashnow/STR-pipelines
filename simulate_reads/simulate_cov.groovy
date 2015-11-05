// Bpipe pipeline to simulate paired end reads from a fasta file

ART='/vlsci/VR0320/hdashnow/art_bin_ChocolateCherryCake'
REF='/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta' 


def get_fname(path) {
    x = path.split("/")[-1]
    return(x)
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
            $ART/art_illumina -i $REF -p -na 
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
// Steps in common
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



/////////////////////////////
// LobSTR

// Paths relative to the analysis directory
PATH_TO_LOBSTR="/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/lobSTR-3.0.3"
LOBSTR_BUNDLE="/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/hg19_v3.0.2"
LOBSTR_REF_PATH="$LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_ref"

load "/vlsci/VR0320/hdashnow/git/STR-pipelines/lobSTR_test/bpipe.config"

//@transform("bam")
// Note the replaceAll with this regex will remove the last extension of the filename

lobSTR = {
    from("fastq.gz","fastq.gz") transform("aligned.bam", "aligned.stats") {
        exec """
            lobSTR \
              --p1 $input1 \
              --p2 $input2 \
              -q --gzip \
              --index-prefix $LOBSTR_REF_PATH/lobSTR_ \
              -o $output.bam.prefix.prefix \
              --rg-sample ${branch.sample} --rg-lib ${branch.lane}
        """, "medium"
    }
}

@filter("sort")
sort_bam = {
    from("aligned.bam") {
        exec "samtools sort $input1 $output.prefix", "medium"
    }
}


@transform("vcf")
allelotype = {
    transform('.bam') to('.lobstr.vcf', '.lobstr.allelotype.stats') {
        def outname = output.prefix[0..-2]
        exec """
            allelotype \
              --command classify \
              --bam $input.bam \
              --noise_model $PATH_TO_LOBSTR/models/illumina_v3.pcrfree \
              --out $output1.prefix \
              --strinfo $LOBSTR_BUNDLE/lobstr_v3.0.2_hg19_strinfo.tab \
              --index-prefix $LOBSTR_REF_PATH/lobSTR_
        """, "medium"
    }
}

/////////////////////////////
// RepeatSeq

// Paths relative to the analysis directory
load "/vlsci/VR0320/hdashnow/git/STR-pipelines/repeatseq/bpipe.config"
// RepeatSeq installation
REPEATSEQ="/vlsci/VR0320/hdashnow/git/repeatseq"

// Reference data
REF="/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta" // Genome
REGIONS="$REPEATSEQ/regions/hg19.2014.regions"
// STRs in the format required by RepeatSeq

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

genotype = {

    doc "Genotype STRs using RepeatSeq"

    from("bam") transform("bam.vcf") {
        exec """
            $REPEATSEQ/repeatseq -emitconfidentsites $input.bam $REF $REGIONS
        """, "medium"
    }
}

/////////////////////////////

// LobSTR
LobSTR_segment = segment {
    lobSTR + sort_bam + index_bam +
    allelotype
}

// RepeatSeq
RepeatSeq_segment = segment {
    align_bwa + index_bam + 
    genotype

}


// Run all

coverages = [5, 10, 15, 20]

run {
    coverages * [
        generate_reads
    ] +
    "%.fq" * [ 
        gzip //+ fastqc
    ] +
    '%_R*.fastq.gz' * [
        set_sample_info +
        align_bwa + index_bam +
        [ allelotype, genotype ]


//        [ LobSTR_segment, RepeatSeq_segment ]
    ] 

}


