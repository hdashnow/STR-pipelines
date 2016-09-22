#!/usr/bin/env Rscript

# Generate a file containing STRs in interesting regions (e.g. coding)
# Generate pathogenic mutants in a single locus or all loci
# Generate normal mutants of all other loci based on the distribution of allele sizes in the reference genome

suppressPackageStartupMessages({
  library('argparse', quietly=TRUE)
})

parser <- ArgumentParser(description="Generate bed file of STR loci to simulate")
parser$add_argument("TRF", metavar="TRFannotation.txt.gz", type="character", nargs=1,
                    help="Gzipped TRF annotation of the genome from UCSC Table Browser")
parser$add_argument("bed_folder", metavar="bed_folder", type="character", nargs=1,
                    help="Location of folder containing annotation of the genome in bed format. Only TRF loci contained within these regions will be used.")
parser$add_argument("-L", "--interval", type="character", default="all",
                    help="Only generate mutants at this single genomic interval, in the form chr1:100-200. [default %(default)s]",
                    metavar="LOCUS")
parser$add_argument("-b", "--background", type="character",
                    help="In the case where -L is used, also mutate all other loci in the normal range and save to this file.",
                    metavar="FILE")
parser$add_argument("-m", "--max", type="integer", default=1000,
                    help="Maximum pathogenic mutant allele to simulate (in repeat units).",
                    metavar="FILE")
parser$add_argument("-O", "--output", type="character", default="stdout",
                    help="Ouput bed file name. [default %(default)s]",
                    metavar="FILE")
# Process command line arguments
args <- parser$parse_args()
#write(str(args), output_bed)

# Input data:
TRF.file = args$TRF #"/Users/hd_vlsci/Documents/reference-data/hg19.simpleRepeat.txt.gz"
ref.folder = args$bed_folder #"/Users/hd_vlsci/Documents/reference-data/str-stats" # containing hg19_RefSeq.*.flat.bed

str.is.integer = function(x) {
  n = as.numeric(x)
  if (is.na(x)) return(FALSE)
  if ( n%%1 == 0) return(TRUE)
  return(FALSE)
}

# Still needs some work. Situation where missing values in any postion can cause errors,
# also start < end not tested for
is.region = function(region) {
  colonsplit = strsplit(region, ':')
  if (length(colonsplit[[1]]) != 2) return(FALSE)
  dashsplit = strsplit(colonsplit[[1]][2], '-')
  if (length(dashsplit[[1]]) != 2) return(FALSE)
  if (str.is.integer(dashsplit[[1]][1]) & str.is.integer(dashsplit[[1]][2])) return(TRUE)
  return(FALSE)
}

if (args$interval == 'all') {
  #output_bed = 'coding_path_STR.bed'
  interval = args$interval
  } else if (!is.region(args$interval)) {
    warning("Region is in the incorrect format. It should be chr1:100-200.")
  } else {
    interval = args$interval
  }

if (args$output == 'stdout') {
  output_bed = stdout()
} else {
  output_bed = args$output
}

max_path = args$max

# Load packages after dealing with commandline input because it takes a while
# (so delay when just using --help)
suppressPackageStartupMessages({
  library('scales', quietly=TRUE)
  library('plyr', quietly=TRUE)
  library('dplyr', quietly=TRUE)
  library('ggplot2', quietly=TRUE)
  library('GenomicRanges', quietly=TRUE)
  library('truncnorm', quietly=TRUE)
})
#ignore.stderr = TRUE # Could be used to get the warnings to print to stderr?

## Functions

parse_multiple_beds = function(files.vect) {
  df = NULL
  for (f in all.files) {
    bed = read.delim(f, header=FALSE,
                     col.names = c("chrom", "start", "end"))
    fname = tail(strsplit(f, '/')[[1]], 1)
    region = substr(fname, 13, nchar(fname) - 9)
    bed$region = region
    df = rbind(df, bed)
  }
  return(df)
}

generate.genotype = function(copyNum, coding.sd) {
  alleles = round(rtruncnorm(2, a=-ceiling(copyNum), b=Inf, mean=0, sd=coding.sd))
  paste(alleles[1], alleles[2], sep='/')
}

generate.uniform.genotypes = function(min.allele, max.allele, ref.allele, n=2) {
  ref.allele = floor(ref.allele)
  alleles = sample(min.allele:max.allele, n, replace=T) - ref.allele
  return(alleles)
}

generate.rand.path.genotype = function(copyNum, coding.sd, max.allele = max_path) {

#  norm.allele = round(rtruncnorm(1, a=-ceiling(copyNum), b=Inf, mean=0, sd=coding.sd))
  norm.allele = ceiling(copyNum) # Set norm.allele to reference allele so don't have to worry about having enough copies of the repeat unit to delete. This is a problem for impure loci.
  path.allele = round(generate.uniform.genotypes(min.allele=copyNum, max.allele,
                                                 ref.allele=copyNum, n=1))
  allele = paste(norm.allele, path.allele, sep='/')
  return(allele)
}

# Write positions and genotypes back to file in bed format
# Assume df has columns: seqnames,start,end,sequence,genotype[,gene]
df.to.bed = function(df,outfile) {
  out.df = df[,c('seqnames','start','end')]
  out.df$start = out.df$start - 1 # To bed format indexing, assuming df has is Granges indexing
  if ('gene' %in% names(df)) {
    out.df$name = paste(df$sequence, df$genotype, df$gene, sep='_')
  } else {
    out.df$name = paste(df$sequence, df$genotype, sep='_')
  }
  write.table(out.df, outfile, sep='\t', quote = F, row.names = F, col.names = F)
}

parse.interval = function(interval) {
  if (!is.region(interval)) return(NA)
  colonsplit = strsplit(interval, ':')
  chrom = colonsplit[[1]][1]
  dashsplit = strsplit(colonsplit[[1]][2], '-')
  return(c(chrom, dashsplit[[1]]))
}


## Main

# Microsatellite data
msats = read.delim(TRF.file, header=FALSE,
                   col.names = c("bin", "chrom", "chromStart", "chromEnd", "name", "period",
                                 "copyNum", "consensusSize", "perMatch", "perIndel", "score",
                                 "A", "C", "G", "T", "entropy", "sequence"))
msats = msats[msats$period<=6,]
msats$length.bp = msats$period * msats$copyNum # Calculate length in bp

# Get RefSeq regions
all.files = list.files(path = ref.folder, pattern = "hg19_RefSeq.*.flat.bed$", full.names = TRUE)
all.refseq = parse_multiple_beds(all.files)

# Get STRs in coding regions (some code from TRF.R)

# BED format is base zero, and the end position is +1 (like python indexing).
# So to convert from bed format to GRanges: start + 1
# sorting seqlevels (i.e. chomosome names) so that all GR objects will have same chromosome ordering.
msats.GR = with(msats, GRanges(chrom, IRanges(chromStart + 1, chromEnd),period = period, copyNum=copyNum, sequence= sequence))
seqlevels(msats.GR) = sort(seqlevels(msats.GR))
msats.GR = sort(msats.GR)
refseq.GR = with(all.refseq, GRanges(chrom, IRanges(start + 1, end), region=region))
seqlevels(refseq.GR) = sort(seqlevels(refseq.GR))
refseq.GR = sort(refseq.GR)

# Annotate region for each STR (that has one)
refseq.GR.by_region = split(refseq.GR, refseq.GR$region)
msats.GR.ann = lapply(refseq.GR.by_region, function(x) {
  # If type is within, the query interval must be wholly contained within the subject interval.
  subsetByOverlaps(msats.GR, x, type = 'within')
})
for (n in names(msats.GR.ann)) {
  #print(n)
  if (length(msats.GR.ann[[n]]) > 0) {
    mcols(msats.GR.ann[[n]])$region = n
  }
}
# Extract out metadata columns
df.list = lapply(msats.GR.ann, mcols)
# Merge together all dfs into a single one
df.all <- ldply(df.list, data.frame)

# Just get the STRs in coding regions
coding.GR = msats.GR.ann[['coding']]

# See what the copy number distribution looks like. Use this for simulation
coding.m = round(median(coding.GR$copyNum))
coding.sd = IQR(coding.GR$copyNum)
coding.df = as.data.frame(coding.GR)


coding.df$genotype = sapply(coding.df$copyNum, generate.genotype, coding.sd)
#df.to.bed(coding.df, 'coding_normal_STR.bed')

### Choose one locus and generate heterzygous pathogenic range genotype
if (interval != 'all') {
  interval.list = parse.interval(interval)
  chrom = interval.list[1]
  start = as.numeric(interval.list[2])
  end = as.numeric(interval.list[3])
  locus.row = which(coding.df$seqnames == chrom & coding.df$start == start & coding.df$end == end)
  path.df = coding.df[locus.row,]
  # Just assuming this will be one row for not, may not generalise
  path.df$genotype = generate.rand.path.genotype(path.df[,'copyNum'], coding.sd, max.allele = max_path)
  df.to.bed(path.df, output_bed)
  # If --background option,
  if (!is.null(args$background)) {
    background.bed = args$background
    background.df = coding.df[-locus.row,]
    df.to.bed(background.df, background.bed)
  }
} else {
  #XXX Could potentially replace the one locus with the pathogenic one, but this isn't in the commandline options currently
  # Replace genotype on this row with a random, potentially pathogenic one
  #coding.df$genotype[locus.row] = generate.rand.path.genotype(coding.df[locus.row,'copyNum'], coding.sd, max.allele = max_path)
  df.to.bed(coding.df, output_bed)
}



#XXX another idea, not coded in options
# ### Generate pathogenic length genotypes at all coding STR loci
#
# # Make a copy of the data frame containing info for all coding STR loci
# random.path = data.frame(coding.df)
# # Generate pathogenic variants of a range of sizes
# random.path$genotype = sapply(random.path$copyNum, generate.rand.path.genotype,
#                               coding.sd, max.allele = max_path)

quit()
