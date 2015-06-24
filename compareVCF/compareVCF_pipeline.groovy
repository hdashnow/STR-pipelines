// Pipeline to compare two vcf files of microsatellite genotypes
// Input: two .vcf files e.g. from different genotypers or samples
// Output: stats on the difference between the two?

// Paths relative to the analysis directory

// Reference data
REF="/vlsci/VR0002/shared/Reference_Files/GATK_bundle_Refs/hg19/ucsc.hg19.fasta"

// Load bpipe.config
load "/vlsci/VR0320/hdashnow/git/STR-pipelines/compareVCF/bpipe.config"

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    }

// -O v for vcf, -O z for vcf.gz
// expecting vcf.gz input?
@filter('trimmed')
trim_variants = {
        exec "bcftools norm -f $REF -O v $input.vcf > $output.vcf"
}

compress_vcf = {
    from('vcf') transform('vcf.gz') {
        exec "bgzip -c $input.vcf > $output.gz"
    }
}

index_vcf = {
    from('vcf.gz') produce(input.gz + '.tbi') {
        exec "tabix -p vcf $input.gz"
    }
    forward input.gz
}

compare_vcfs = {
    from('vcf.gz', 'vcf.gz') produce(input1.vcf.prefix + '-' + input2.vcf.prefix + '.diff.txt') { 
        exec """
            vcf-compare -g $input1.gz $input2.gz > $output.txt
        """
    }
}

run {
    "%.vcf" * [
        set_sample_info + trim_variants +
        compress_vcf + index_vcf
    ] + compare_vcfs
}
