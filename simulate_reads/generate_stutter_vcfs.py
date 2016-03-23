#!/usr/bin/env python
"""Generates VCFs to simulate STR variants and stutter.
To be fed into genome and read simulation stages.

Output:
A single vcf file containing the true genotype that is being simulated.
A series of vcf files for the stutter (one for each STR variant)
A bed file corresponding to the region around the STR.
A single file with the list of all vcf files and their desired frequencies
"""

import sys
from argparse import (ArgumentParser, FileType)
from Bio.Seq import Seq #BioPython
from Bio.Alphabet import generic_dna
import vcf #PyVCF

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = ArgumentParser(description='Produce VCFs and corresponding stutter frequencies for a given set of STR loci')
    parser.add_argument(
        '--ref', type=str, required=True,
        help='Fasta reference')
    parser.add_argument(
        '--bed', type=str, required=True,
        help='bed file containing genomic locations of STRs and their repeat units. Genomic locations should be relative to the fasta reference.')
    parser.add_argument(
        '--output', type=str, required=False,
        help='Base name for output files')
    parser.add_argument(
        '--base0', action='store_true',
        help='Genomic positions in bed file and region are 0-based. Otherwise assumed to be 1-based.')
    parser.add_argument(
        '--flank', type=int, default=100000,
        help='Number of flanking base to include in the output bed file on either side of the STR.')
    return parser.parse_args()

def circular_permuted(x):
    """Generates all possible circular permutations of the input.

    Args:
        x (str or any iterable?)

    Returns:
        list: All circular permutations of x
    """
    return([x[i:] + x[:i] for i in range(len(x))])

def self_and_rev_complement(in_dna):
    """Returns the input DNA sequence and its reverse complement

    Args:
        in_dna (str): valid DNA bases e.g. ACTGN

    Returns:
        list: [in_dna, reverse_complement]
    """
    all_possible = [in_dna]
    # Get reverse complement
    dna = Seq(in_dna, generic_dna)
    rev_complement = str(dna.reverse_complement())
    all_possible.append(rev_complement)
    return(all_possible)

def normalise_str(in_dna):
    """Find all possible eqivalent Short Tandem Repeat (STR) DNA sequences.
    And return the first alphabetically.
    For example, TA = AT. But would return AT.

    Args:
        in_dna (str): STR repeat unit consisting of valid DNA bases e.g. ACTGN

    Returns:
        str: normalised STR sequence
    """
    all_possible = []
    # Circularly permute original sequence and reverse complement
    for seq in self_and_rev_complement(in_dna):
        for permuted_seq in circular_permuted(seq): # Switch to faster permutation (6)
            all_possible.append(permuted_seq)

    # Sort and take the first
    all_possible.sort()
    return(all_possible[0])

def parse_bed(bedfilename, bed_dict = {}, position_base = 1):
    """Parse regions from bed file. Ignore lines starting with #.

    Args:
        bedfilename (str): bed format text file in in format: chr start stop [name] ... (additional columns ignored)
        bed_dict (dict): existing genomic regions to include in the output.
            format: bed_dict[unique_id] = {"chr":str, "start":int, "stop":int, "name":None}
        position_base (int): 0 or 1. The starting position the genomic regions
            are measured relative to. If 1, all postions will be converted to 0.

    Returns:
        dict: bed_dict[unique_id] = {"chr":str, "start":int, "stop":int, "name":None}
            Genomic regions, in base-0.
    """
    with open(bedfilename) as bedfile:
        for bedfile_line in bedfile:
            if bedfile_line.startswith("#"):
                continue
            split_line = bedfile_line.split()
            ref_chr = split_line[0]
            ref_start = int(split_line[1]) - position_base
            ref_stop = int(split_line[2]) - position_base
            # If repeat start and end are the wrong way around, swap them
            if ref_stop < ref_start:
                ref_start = int(split_line[2]) - position_base
                ref_stop = int(split_line[1]) - position_base
                sys.stderr.write('Warning, bed start position greater than end position for line:')
                sys.stderr.write(bedfile_line)
            if len(split_line) > 3:
                name = split_line[3] #XXX Need to parse out STR repeat unit here
            else:
                name = None
            unique_id = "{0}:{1}-{2}".format(ref_chr,ref_start,ref_stop)
            if unique_id not in bed_dict:
                bed_dict[unique_id] = {"chr":ref_chr, "start":ref_start, "stop":ref_stop, "name":name}
    return bed_dict

def get_vcf_writer(vcf_outfile):
    """Generate a vcf writer object.
    Writes a template vcf containing header info (can be deleted afterwards?)"""

    template = 'template.vcf'
    with open(template, "w") as vcf_template:
        vcf_template.write('##fileformat=VCFv4.1\n')
        vcf_template.write('##source={}\n'.format('generate_stutter_vcfs.py'))
        vcf_template.write('##reference={}\n'.format('ucsc.hg19.fasta'))
        vcf_template.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat Unit">\n')
        vcf_template.write('##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference Length of Repeat">\n')
        vcf_template.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf_template.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')

    with open(template) as vcf_template:
        vcf_reader = vcf.Reader(vcf_template)
        vcf_writer = vcf.Writer(open(vcf_outfile, 'w'), vcf_reader)
    return vcf_writer

def main():
    # Parse command line arguments
    args = parse_args()
    outfile_base = args.output

    vcf_out = outfile_base + '.vcf'

    if args.base0:
        position_base = 0
    else:
        position_base = 1

    bed_dict = {}
    #bed_dict = parse_bed(args.bed, bed_dict, position_base)

    vcf_writer = get_vcf_writer(vcf_out)

    record = vcf.model._Record(CHROM='chr6', POS=16243709, ID='.', REF='T',
                                ALT=[vcf.model._Substitution('TGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT')],
                                QUAL='.', FILTER='PASS', INFO={'RU':'GCT', 'RL':29},
                                FORMAT='.', sample_indexes=[], samples=None)
    vcf_writer.write_record(record)

if __name__ == '__main__':
    main()