// bpipe pipeline to run STR_FM
// Steps:
// - Detect location of all STR in reference genome
// - Profile STR from short reads using STR_FM

////////////////////////////////////////////////////
// Detect location of all STR in reference genome //
////////////////////////////////////////////////////


// Requirement: reference genome in FASTA format --> REF
BASE="/mnt/storage/harrietd/git"
REF="/mnt/storage/shared/genomes/hg19/fasta/hg19.fa"
REFNAME="hg19"
STR_FM="$BASE/STR-FM"
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

run {
    motif_lens * [ detect_ref_STRs + format_STRs ]
}

