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

compare_vcfs = {
    from('vcf', 'vcf') produce(input1.vcf.prefix + '-' + input2.vcf.prefix + '.diff.indv_in_files', 
        input1.vcf.prefix + '-' + input2.vcf.prefix + '.diff.sites_in_files', 
        input1.vcf.prefix + '-' + input2.vcf.prefix + '.summary.txt') {
        exec """
            vcftools --vcf $input1 --diff $input2 --diff-site-discordance
            --out $output1.prefix.prefix 2> $output.txt
        """
    }
}

run {
    "%.vcf" * [
        set_sample_info + trim_variants 
    ] + compare_vcfs
}
