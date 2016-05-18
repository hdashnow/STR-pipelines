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
import pandas as pd

__author__ = 'Harriet Dashnow'
__credits__ = ['Harriet Dashnow']
__license__ = 'MIT'
__version__ = '0.1.0'
__email__ = 'h.dashnow@gmail.com'

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = ArgumentParser(description=('Produce VCFs and corresponding stutter'
    ' frequencies for a given set of STR loci. Also provides a bed file for each'
    ' locus defining a region around that locus.'))
    parser.add_argument(
        'ref', type=str, #XXX switch to positional argument?
        help='Fasta reference')
    parser.add_argument(
        'bed', type=str, #XXX switch to positional argument?
        help='bed file containing genomic locations of STRs and their repeat units. Genomic locations should be relative to the fasta reference.')
    parser.add_argument(
        '--output', type=str, required=False,
        help='Base name for output files, including vcfs and bed files.') #XXX should have a default for this?
    parser.add_argument(
        '--truth', type=str, required=False,
        help='File name for output vcf of true genotypes for all loci. (default: output base name .truth.vcf)')
    parser.add_argument(
        '--stutter_output', type=str, required=False,
        help='File giving names of stutter vcf files with corresponding stutter probabilities. (default: stdout)')
    parser.add_argument(
        '--base0', action='store_true',
        help='Genomic positions in bed file and region are 0-based. Otherwise assumed to be 1-based.')
    parser.add_argument(
        '--flank', type=int, default=10000,
        help='Number of flanking base to include in the output bed file on either side of the STR. (default: %(default)s)')
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
                    repeatunit = name.upper()
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
                ref_seg = bases_to_right[j:j+repeatunitlen]
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

def combine_stutter(deltas1, probs1, deltas2, probs2, rescale_probs = True):
    """Combine the stutter distributions from two different alleles to give the stutter distrobutions for an entire genotype.

    Args:
        deltas1, deltas2 (list):
        probs1, probs2 (list):
        rescale_probs (bool): If the individual input probs don't add up to 1, rescale them so they do before combining them.

    Returns:
        (deltas, probs)
        deltas (list): Combinded deltas
        probs (list): Combined probs. Rescaled to sum to 1.
    """
    allele1 = pd.DataFrame(data=probs1, index=deltas1, columns=["probs1"])
    allele2 = pd.DataFrame(data=probs2, index=deltas2, columns=["probs2"])

    both = allele1.join(allele2, how='outer').fillna(0)

    if rescale_probs: # rescale individual allele probs to add to 1
        both["probs1"] = both["probs1"]/sum(both["probs1"])
        both["probs2"] = both["probs2"]/sum(both["probs2"])

    combined = both["probs1"] + both["probs2"]
    combined = combined/sum(combined) # rescale to add to 1

    return list(combined.index), list(combined)

def main():
    # Parse command line arguments
    args = parse_args()
    outfile_base = args.output

    if args.truth:
        truth_fname = args.truth
    else:
        truth_fname = outfile_base + '.truth.vcf'
    vcf_truth = get_vcf_writer(truth_fname)

    if args.stutter_output:
        vcf_probs_writer = open(args.stutter_output, "w")
    else:
        vcf_probs_writer = sys.stdout

    if args.base0:
        position_base = 0
    else:
        position_base = 1

    if args.seed:
        random.seed(args.seed)

    # Hard-coding some probabilities for testing
    # Gymrek, M. (2016). PCR-free library preparation greatly reduces stutter noise at short tandem repeats.
    # Retrieved from http://biorxiv.org/lookup/doi/10.1101/043448
    #XXX But I'm not sure about these numbers. They don't add up to 1 for 3 and 4 bp repeat units - not being used currently.
    # probability of stutter occuring:
    prop_of_stutter = {1: 0.17, 2: 0.038, 3: 0.011, 4: 0.0069, 5: 0.0074, 6: 0.014}
    # if stutter occurs, distribution of deltas:
    RU2_stutter = [[-3, -2, -1, 1, 2], [0.02, 0.2, 0.6, 0.09, 0.02]]
    RU3_stutter = [[-3, -2, -1, 1, 2], [0.01, 0.1, 0.45, 0.05, 0.01]]
    RU4_stutter = [[-2, -1, 1, 2], [0.05, 0.34, 0.04, 0.005]]

    # Parse STR regions that need to be simulated
    bed_dict = parse_bed(args.bed, position_base)
    # get corresponding bit of the fasta file
    fastafile = pysam.Fastafile(args.ref)
    for region in bed_dict:
        chrom = bed_dict[region]['chr']
        start = bed_dict[region]['start']
        stop = bed_dict[region]['stop']
        name = bed_dict[region]['name']
        repeatunit = bed_dict[region]['repeatunit']
        ref_sequence = fastafile.fetch(chrom, start, stop).upper()  # note zero-based indexing c.f. 1-based indexing in vcf
        # XXX need to switch back to 1-based indexing when writing to the vcf! Make a function for this?

        # Write a bed file of the bases around the given region
        #XXX Note, this is writing out the bed file in base 0, change to base 1 (or give option to)? Give info in help about this?
        bed_out = '{}_{}.bed'.format(outfile_base,region)
        with open(bed_out, "w") as f:
            f.write('{0}\t{1}\t{2}\t{3}\n'.format(chrom, start - args.flank, stop + args.flank, name))

        #delta = 1 #XXX need to generate delta, or get from input

        # Generate a genotype for these - totally random, or heterozygous pathogenic?
        allele1_delta = 10
        allele2_delta = 0
        allele1 = mutate_str(ref_sequence, repeatunit, delta = allele1_delta)
        allele2 = mutate_str(ref_sequence, repeatunit, delta = allele2_delta)

        # Write the true alleles (the basis for the stutter simulation)
        # XXX need to figure out how to write genotype, not just two alleles?
        # XXX also shouldn't include ALT for allele if same as ref - should be in genotype instead
        record = vcf.model._Record(CHROM=chrom, POS=start, ID='.', REF=ref_sequence,
                    ALT=[vcf.model._Substitution(allele1), vcf.model._Substitution(allele2)],
                    QUAL='.', FILTER='PASS', INFO={'RU':repeatunit},
                    FORMAT='.', sample_indexes=[], samples=None)
        vcf_truth.write_record(record)

        # Calculate stutter probability profile for each allele
        # Parameters: repeat unit size, repeat length?
        deltas1 = [-3, -2, -1, 0, 1, 2]
        probs1 =  [0.2, 0.3, 0.5, 0.7, 0.4, 0.2]
        deltas2 = [7, 8, 9, 10, 11, 12]
        probs2 =  [0.2, 0.3, 0.5, 0.7, 0.4, 0.2]
        stutter_deltas,stutter_probs = combine_stutter(deltas1, probs1, deltas2, probs2, rescale_probs = True)

        # Generate stutter alleles
        for delta, prob in zip(stutter_deltas,stutter_probs):
            stutter_fname = outfile_base + '.stutter_{0}.vcf'.format(delta)
            vcf_stutter = get_vcf_writer(stutter_fname)
            mutatant_allele = mutate_str(ref_sequence, repeatunit, delta = delta)
            if delta == 0:
                record = vcf.model._Record(CHROM=chrom, POS=start, ID='.', REF=ref_sequence,
                            ALT=[],
                            QUAL='.', FILTER='PASS', INFO={'RU':repeatunit},
                            FORMAT='.', sample_indexes=[], samples=None)
            else:
                record = vcf.model._Record(CHROM=chrom, POS=start, ID='.', REF=ref_sequence,
                            ALT=[vcf.model._Substitution(mutatant_allele)],
                            QUAL='.', FILTER='PASS', INFO={'RU':repeatunit},
                            FORMAT='.', sample_indexes=[], samples=None)
            vcf_stutter.write_record(record)
            # write the filename and corresponding stutter probability for use in later pipeline stages
            vcf_probs_writer.write('{0}\t{1}\t{2}\n'.format(stutter_fname, prob, bed_out))

        # Write true genotype vcf, stutter vcfs and their corresponding probabilities



    #vcf_writer.write_record(record)




if __name__ == '__main__':
    main()
