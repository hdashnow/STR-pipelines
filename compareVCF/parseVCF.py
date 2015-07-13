# parse VCF files of STR genotype calls (from LobSTR and RepeatSeq)
# Compare the results

import vcf
import pandas as pd
import sys

dir = '/Users/hd_vlsci/Documents/git/STR-pipelines/data/intersections_LobSTR-RepeatSeq/'
lobstr = dir + 'intersection0_10.vcf' # Called by LobSTR only
both = dir + 'intersection0_10_1.vcf' # Called by LobSTR and RepeatSeq
repeatseq = dir + 'intersection0_11.vcf' # Called by RepeatSeq only

vcf_reader = vcf.Reader(open(repeatseq, 'r'))

sample1 = vcf_reader.samples[0] # Get list of samples in the VCF. save the first

for record in vcf_reader:
    print(pd.DataFrame(record.INFO))
    print(record.genotype(sample1)['GT'])
    sys.exit()


# Data to extract:
# Chrom
# Pos
# Genotype
GT
# Allele Length Offset(s)" - need to figure out what this is!
# Depth
DP
# Repeat unit
RU
# Reference length of repeat
RL

# Things to calculate:
