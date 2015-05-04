// bpipe pipeline to run STR_FM
// Steps:
// - Detect location of all STR in reference genome
// - Profile STR from short reads using STR_FM

///////////////////////////
// Requirements/settings //
///////////////////////////

// reference genome in FASTA format --> REF
BASE="/mnt/storage/harrietd/git"
REF="/mnt/storage/shared/genomes/hg19/fasta/hg19.fa"
REFNAME="hg19"

// Path to STR-FM scripts
STR_FM="$BASE/STR-FM"

// index of mapping program (bwa, bowtie, etc) 
// reference genome in FASTA and in 2bit file --> /path/to/2bit/ref.2bit (use utility from UCSC genome browser to create 2bit file version of reference genome)
REF2bit = "2bitref.2bit"

// STR error rates (can be downloaded from https://usegalaxy.org/u/guru%40psu.edu/h/error-rates-files) --> errorrate.bymajorallele

GALAXY_TOOLS="/mnt/storage/harrietd/git/galaxy/tools"

////////////////////////////////////////////////////
// Detect location of all STR in reference genome //
////////////////////////////////////////////////////

motif_lens=["mono", "di", "tri", "tetra"]

detect_ref_STRs = {
    doc "detect STR in reference genome"
    output.dir = "detectSTRs"

    branch.outname = branch.name

    if (branch.name == "mono") {
        branch.period = 1
        branch.minlength = 4
    }
    if (branch.name == "di") {
        branch.period = 2
        branch.minlength = 6
    }
    if (branch.name == "tri") {
        branch.period = 3
        branch.minlength = 6
    }
    if (branch.name == "tetra") {
        branch.period = 4
        branch.minlength = 8
    }

    produce("${REFNAME}.${outname}.out") {

        exec """
            python $STR_FM/microsatellite.py $REF --fasta --period=$period --partialmotifs --minlength=$minlength --prefix=0 --suffix=0 --hamming=0 --multipleruns --flankdisplay=0  --splitbyvalidity  > $output
        """
    }
}

format_STRs = {
    doc "format detected STRs"
    output.dir = "detectSTRs"

    exec """
        cat $input.out | awk 'BEGIN{FS="\t";OFS="\t"};{print \$6,\$2,\$2+\$1,\$4,\$1,length(\$4) }' > $output.TR
    """
}

///////////////////////////////////////////////
// Profile STR from short reads using STR-FM //
///////////////////////////////////////////////

// Requirements:
// fastq input in sangerfq Phred scale --> ${INPUT}.fastq
// location of all STR in reference genome (use PBS script name "sampleSTR_reference_profiling.txt) --> /path/to/STR/in/reference/genome.TR (you can make 4 separated TR files for 4 types of STRs)

// ref=/path/to/reference/sequence/and/bwa/index/ref.fa
// export PYTHONPATH=/path/to/galaxy-dist/lib/
// galaxydir=/path/to/galaxy-dist/tools

detect_read_STR = {
    doc "detect STR in short read"
    output.dir = "detect_reads"

    produce("${sample}.${outname}.out") {
        exec """
             gzip -cd $input.fastq.gz | python $STR_FM/microsatellite.py  --fastq --period=$period --partialmotifs --minlength=5 --prefix=20 --suffix=20 --hamming=0 --multipleruns  > $output
        """
    }
}

@transform("new")
rename_reads = {
    doc "change read name (replace space with underscore)"
    output.dir = "detect_reads"

    exec """
      python $STR_FM/changespacetounderscore_readname.py $input.out  $output.new 6
    """
}

flanking_seq = {
    doc "fetch flanking sequence of STRs for mapping"
    output.dir "detect_reads"

    produce("${output}_ff_L.txt", "${output}_ff_R.txt"){
        exec """
            python $STR_FM/pair_fetch_DNA_ff.py $input.new $outputs 20 20
        """
    }
}

align_bwa = { 

    doc "map flanking sequences using BWA - uniquely mapped no indel no deletion"
    output.dir "align"

    // transform two .txt fils (L and R flanks) into two .sai files and a .bam file
    from("txt","txt") transform("sai","sai","bam") {

        // Step 1 - run both bwa aln commands in parallel
        multi "bwa aln -n 0 -o 0 $REF $input1 > $output1",
              "bwa aln -n 0 -o 0 $REF $input2 > $output2"

        // Step 2 - bwa sampe
        exec """
            bwa sampe $REF $output1 $output2 $input1 $input2 |
            samtools view -Sbu -F 12 -q 1 - > $output.bam
        """
    }
}

sort_bam = {
    doc "sort result by read name and fix header"
    output.dir "align"

    produce("${input.prefix}.sorted.sam") {
        exec """
            samtools sort -on $input.bam |
            samtools view -h -o $output.sam -
         """
    }
}

@transform("RF")
merge_flanks = {
    doc "merge faux paired end reads"
    output.dir "align"

    exec """
        python $STR_FM/PEsortedSAM2readprofile.py $input.sam $REF2bit 100 250 $output.RF
     """
}


// Not sure how to actually use this join function as the instructions has it joining all the files together using an x = x + 1 style strategy
join_mapped = {
    doc "join mapped coordinate with STR length using read name"
    output.dir "align"

    from("new","RF") transform("j") {
        exec """
            python $GALAXY_TOOLS/filters/join.py $input.new $input.RF 6 1 $output.j "" "" --index_depth=3 --buffer=50000000 --fill_options_file='None'
         """
    }
}


run {
    motif_lens * [ 
        detect_ref_STRs + format_STRs + 
        "%.fastq.gz" * [ 
            detect_read_STR + rename_reads + flanking_seq + align_bwa + sort_bam
        ]
    ]
}

