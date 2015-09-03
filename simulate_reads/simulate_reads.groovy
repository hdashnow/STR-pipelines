// Bpipe pipeline to simulate paired end reads from a fasta file

ART='/vlsci/VR0320/hdashnow/art_bin_ChocolateCherryCake'
//REF='/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta' 

generate_reads = {
    def samplename = branch.name
    produce(samplename + '_L001_R1.fq', samplename + '_L001_R2.fq') {
        def outname = output.prefix[0..-2]
        exec """
            $ART/art_illumina -i $input.fasta -p -na 
                -l 150 -ss HS25 -f 20 
                -m 500 -s 50 
                -o $outname
        """
    }
}

gzip = {
    transform('.fq') to('.fq.gz') {
        exec "gzip -c $input.fq > $output.gz" 
    }
}

run {
    "%.fasta" * [
        generate_reads
    ] +
    "%.fq" * [ 
        gzip 
    ]
}
