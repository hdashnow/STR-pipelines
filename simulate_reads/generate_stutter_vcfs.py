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
import pysam #note: can be tricky to install, used bioconda channel
import random

__author__ = 'Harriet Dashnow'
__credits__ = ['Harriet Dashnow']
__license__ = 'MIT'
__version__ = '0.1.0'
__email__ = 'h.dashnow@gmail.com'

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
    parser.add_argument(
        '--seed', required=False,
        help='Random seed (can be any hashible input).')
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

def is_dna(a):
    """Return True if input contains only DNA characters, otherwise return False

    Args:
        a (str): string to test

    Returns:
        bool: True if all characters are DNA, otherwise False
    """
    dna_chars = 'atcgnATCGN'
    return all(i in dna_chars for i in a)

def parse_bed(bedfilename, position_base = 1, bed_dict = {}):
    """Parse regions from bed file. Ignore lines starting with #.

    Args:
        bedfilename (str): bed format text file in in format: chr start stop [name] ... (additional columns ignored)
        bed_dict (dict): existing genomic regions to include in the output.
            format: bed_dict[unique_id] = {'chr':str, 'start':int, 'stop':int, 'name':None}
        position_base (int): 0 or 1. The starting position the genomic regions
            are measured relative to. If 1, all postions will be converted to 0.

    Returns:
        dict: bed_dict[unique_id] = {'chr':str, 'start':int, 'stop':int, 'name':None}
            Genomic regions, in base-0.
    """
    with open(bedfilename) as bedfile:
        for bedfile_line in bedfile:
            if bedfile_line.startswith('#'):
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
                if is_dna(name):
                    repeatunit = name
            else:
                name = None
                repeatunit = None
            unique_id = '{0}:{1}-{2}'.format(ref_chr,ref_start,ref_stop)
            if unique_id not in bed_dict:
                bed_dict[unique_id] = {'chr':ref_chr, 'start':ref_start,
                                        'stop':ref_stop, 'name':name,
                                        'repeatunit':repeatunit}
    return bed_dict

def mutate_str(ref_sequence, repeatunit, delta):
    """Mutate a DNA sequence containing a microsatellite by inserting or
    deleating repeat units of that microsatellite.

    Args:
        ref_sequence (str): The reference DNA sequence to be mutated.
        repeatunit (str): DNA repeat unit to be inserted/deleted from the
            ref_sequence.
        delta (int): Number of repeat units to add/remove. Positive if creating
            an insertion, negative if a deletion. 0 will return the input
            sequence.

    Returns:
        str: The mutated DNA sequence.
    """
    # Check ref_sequence and repeatunit are DNA
    if not is_dna(ref_sequence):
        raise ValueError("{0} is not a valid DNA sequence".format(ref_sequence))
    if not is_dna(repeatunit):
        raise ValueError("{0} is not a valid DNA sequence".format(repeatunit))
    repeatunitlen = len(repeatunit)
    if delta == 0:
        return(ref_sequence)
    if delta > 0: # Insertion
        # select a position at random, sampling without replacement
        i_max = len(ref_sequence)-repeatunitlen
        for i in random.sample(range(i_max), i_max):
            # Check that there is a repeat unit to the right of this position
            bases_to_right = ref_sequence[i:i+repeatunitlen]
            if normalise_str(bases_to_right) == normalise_str(repeatunit):
                # Insert repeat units at that position
                # (by replicating the reference version of the repeat unit)
                new_sequence = ref_sequence[:i] + bases_to_right * delta + ref_sequence[i:]
                return(new_sequence)
        # Check the last few bases for the repeat unit
        last_bases = ref_sequence[-repeatunitlen:]
        if normalise_str(last_bases) == normalise_str(repeatunit):
            # Insert repeat units at end of sequence
            new_sequence = ref_sequence + last_bases * delta
            return(new_sequence)
        # Check to make sure a solution was found
        raise ValueError("The repeat unit {0} was not found in {1}".format(repeatunit, ref_sequence))
    if delta < 0: # Deletion
        deletion_size = repeatunitlen * -delta
        if deletion_size > len(ref_sequence):
            raise ValueError("Deletion of {0} {1} repeat units is larger than the input sequence {2}".format(-delta, repeatunit, ref_sequence))
        i_max = len(ref_sequence) - deletion_size + 1
        for i in random.sample(range(i_max), i_max):
            bases_to_right = ref_sequence[i:i+deletion_size]
            # Check if the first few bases contain the repeat unit (or a transposition of it)
            for j in range(len(bases_to_right)):
                ref_seg = bases_to_right[j:j+2]
                # Find the version of the repeat unit present in the ref
                if normalise_str(repeatunit) == normalise_str(ref_seg):
                    ref_unit = ref_seg
                    # Check if there are enough repeat units to be deleted
                    if bases_to_right == ref_unit * -delta:
                        # Generate mutated sequence
                        new_sequence = ref_sequence[:i] + ref_sequence[i+deletion_size:]
                        return(new_sequence)
                else:
                    break
        raise ValueError("There were not {0} copies of {1} repeat unit available to be deleted in {2}.".format(-delta, repeatunit, ref_sequence))

def get_vcf_writer(vcf_outfile):
    """Generate a vcf writer object.
    Writes a template vcf containing header info (can be deleted afterwards?)"""

    template = 'template.vcf'
    with open(template, 'w') as vcf_template:
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

    if args.seed:
        random.seed(args.seed)

    # Parse STR regions that need to be simulated
    bed_dict = parse_bed(args.bed, position_base)
    print(bed_dict)
    # get corresponding bit of the fasta file
    fastafile = pysam.Fastafile(args.ref)
    for region in bed_dict:
        chrom = bed_dict[region]['chr']
        start = bed_dict[region]['start']
        stop = bed_dict[region]['stop']
        name = bed_dict[region]['name']
        repeatunit = bed_dict[region]['repeatunit']
        ref_sequence = fastafile.fetch(chrom, start, stop)  # note zero-based indexing c.f. 1-based indexing in vcf

    # Generate a genotype for these - totally random, or heterozygous pathogenic?
    #variant_str = mutate_str(ref_sequence, repeatunit, delta)

    # Calculate stutter probability profile for each allele
    # Parameters: repeat unit size, repeat length?

    # Generate stutter for each allele

    # Write true genotype vcf, stutter vcfs and their corresponding probabilities

    vcf_writer = get_vcf_writer(vcf_out)

    record = vcf.model._Record(CHROM='chr6', POS=16243709, ID='.', REF='T',
                                ALT=[vcf.model._Substitution('TGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT')],
                                QUAL='.', FILTER='PASS', INFO={'RU':'GCT', 'RL':29},
                                FORMAT='.', sample_indexes=[], samples=None)
    #vcf_writer.write_record(record)




if __name__ == '__main__':
    main()
