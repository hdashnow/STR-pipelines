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
        'ref', type=str,
        help='Fasta reference')
    parser.add_argument(
        'bed', type=str,
        help='bed file containing genomic locations of STRs and their repeat units. Genomic locations should be relative to the fasta reference.')
    parser.add_argument(
        '--stutter', type=str, required=False,
        help='csv file containing stutter models.') #XXX need to describe this properly
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
        bedfilename (str): bed format text file in in
        format: chr start stop name [...] (additional columns ignored)
        where name is in the format repeatunit_genotype, e.g.
        CAG_-2/1 (CAG repeat unit, 2 CAG deleation/1 CAG insertion)
        AT_0/3 (AT repeat unit, same as reference/3 AT insertion)
        bed_dict (dict): existing genomic regions to include in the output.
            format: bed_dict[unique_id] = {'chr':str, 'start':int, 'stop':int, 'name':None}
        position_base (int): 0 or 1. The starting position the genomic regions
            are measured relative to in the input file. If 0, all postions will
            be converted to 1. Assumed to be 1 by default.

    Returns:
        dict: bed_dict[unique_id] = {'chr':str, 'start':int, 'stop':int,
                                    'name':str or None, 'repeatunit':str or None,
                                    'deltas': [int, int] or None}
            Genomic regions, in base-1.
    """
    with open(bedfilename) as bedfile:
        if position_base == 0:
            base_shift = 1
        else:
            base_shift = 0
        for bedfile_line in bedfile:
            if bedfile_line.startswith('#') or bedfile_line == '\n':
                continue
            split_line = bedfile_line.split()
            ref_chr = split_line[0]
            ref_start = int(split_line[1]) + base_shift
            ref_stop = int(split_line[2]) + base_shift
            # If repeat start and end are the wrong way around, swap them
            if ref_stop < ref_start:
                ref_start = int(split_line[2]) - position_base
                ref_stop = int(split_line[1]) - position_base
                sys.stderr.write('Warning, bed start position greater than end position for line:')
                sys.stderr.write(bedfile_line)
            if len(split_line) > 3:
                name = split_line[3] #XXX Need to parse out STR repeat unit here
                split_name = name.split('_')
                if len(split_name) >= 2:
                    repeatunit = split_name[0]
                    if is_dna(repeatunit):
                        repeatunit = repeatunit.upper()
                    genotype = split_name[1].split('/')
                    if len(genotype) == 2:
                        deltas = [int(x) for x in genotype]
            else:
                name = None
                repeatunit = None
                deltas = None
            unique_id = '{0}-{1}-{2}'.format(ref_chr,ref_start,ref_stop)
            if unique_id not in bed_dict:
                bed_dict[unique_id] = {'chr':ref_chr, 'start':ref_start,
                                        'stop':ref_stop, 'name':name,
                                        'repeatunit':repeatunit,
                                        'deltas': deltas}
    return bed_dict

def parse_stutter(filename):
    """Parse a csv file of microsatellite stutter frequencies. The first column
    should contain the repeat unit length (no header needed), the remaining
    columns should be named according to the deltas (number of repeat units
    added/deleted e.g. -1 means 1 repeat unit deleated). The values are the
    frequencies of each of these stutter variants.

    Args:
        filename (str): csv file

    Returns:
        pandas.DataFrame: Column names are deltas, rows are repeat unit lengths.
    """
    with open(filename, 'rt') as f:
        data_table = pd.read_csv(f, index_col=0)
    return(data_table)

def mutate_str(ref_sequence, repeatunit, delta, random=False):
    """Mutate a DNA sequence containing a microsatellite by inserting or
    deleating repeat units of that microsatellite.

    Args:
        ref_sequence (str): The reference DNA sequence to be mutated.
        repeatunit (str): DNA repeat unit to be inserted/deleted from the
            ref_sequence.
        delta (int): Number of repeat units to add/remove. Positive if creating
            an insertion, negative if a deletion. 0 will return the input
            sequence.
        random (bool): Insert/delete repeat units from a random position.
        If False, insert in the left-most position. Note, GATK seems to silently
        fail on vcfs with indels to the right of an imperfect repeat.

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
        if random:
            i_range = random.sample(range(i_max), i_max)
        else:
            i_range = range(i_max)
        for i in i_range:
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
        if random:
            i_range = random.sample(range(i_max), i_max)
        else:
            i_range = range(i_max)
        for i in i_range:
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

def trim_indel(ref, alt):
    """Generate the shortest representation of an indel my removing the rightmost bases.
    XXX Currently only working with a single alt, should be generalised for multiple.

    Args:
        ref (str): The version of the sequence in the reference genome.
        alt (str): Alternative version of the sequence caused by an insertion or
        deletion.

    Returns:
        (ref_normalised str, alt_normalised str)
    """
    if min(len(ref), len(alt)) <= 1:
        return ref, alt
    for i in range(1, min(len(ref), len(alt))):
        if ref[-i] != alt[-i]: # Check if the rightmost bases is different
            return ref[:-i+1], alt[:-i+1]
    return ref[:-i], alt[:-i]

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

    if args.stutter:
        stutterDF = parse_stutter(args.stutter)
    else:
        stutterDF = parse_stutter('stutter_model.csv') #XXX This requires that I know where this file is. so maybe refer to it relative to the location of this script?

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

    # Parse STR regions that need to be simulated
    bed_dict = parse_bed(args.bed, position_base)
    # get corresponding bit of the fasta file
    fastafile = pysam.Fastafile(args.ref)
    for region in bed_dict:
        chrom = bed_dict[region]['chr']
        start = bed_dict[region]['start'] # These positions are in base-1
        stop = bed_dict[region]['stop']
        name = bed_dict[region]['name']
        repeatunit = bed_dict[region]['repeatunit']
        deltas = bed_dict[region]['deltas']
        # note fetch step requires zero-based indexing c.f. 1-based indexing in vcf and bed dict
        ref_sequence = fastafile.fetch(chrom, start - 1, stop - 1).upper()

        # Write a bed file of the bases around the given region
        bed_out = '{}.{}.bed'.format(outfile_base,region)
        with open(bed_out, "w") as f:
            f.write('{0}\t{1}\t{2}\t{3}\n'.format(chrom, start - args.flank, stop + args.flank, name))

        #delta = 1 #XXX need to generate delta, or get from input

        # Generate a genotype for these - totally random, or heterozygous pathogenic?
        if deltas:
            allele1_delta = deltas[0]
            allele2_delta = deltas[1]
        else:
            # XXX make this some kind of random number generator
            allele1_delta = 10
            allele2_delta = 0
        allele1 = mutate_str(ref_sequence, repeatunit, delta = allele1_delta)
        allele2 = mutate_str(ref_sequence, repeatunit, delta = allele2_delta)

        # Get stutter from stutterDF based on the repeatunit length
        probs1 = stutterDF.iloc[len(repeatunit)].values
        probs2 = stutterDF.iloc[len(repeatunit)].values
        deltas1 = stutterDF.columns.values.astype(int) + allele1_delta
        deltas2 = stutterDF.columns.values.astype(int) + allele2_delta
        stutter_deltas,stutter_probs = combine_stutter(deltas1, probs1, deltas2, probs2, rescale_probs = True)
        # Possible extension:
        # Calculate stutter probability profile for each allele using a forumlar instead of observed data
        # Parameters: repeat unit size, repeat length?

        # Write the true alleles (the basis for the stutter simulation)
        # XXX need to figure out how to write genotype, not just two alleles?
        # XXX also shouldn't include ALT for allele if same as ref - should be in genotype instead
        record = vcf.model._Record(CHROM=chrom, POS=start, ID='.', REF=ref_sequence,
                    ALT=[vcf.model._Substitution(allele1), vcf.model._Substitution(allele2)],
                    QUAL='.', FILTER='PASS', INFO={'RU':repeatunit},
                    FORMAT='.', sample_indexes=[], samples=None)
        vcf_truth.write_record(record)

        # Generate stutter alleles
        for delta, prob in zip(stutter_deltas,stutter_probs):
            stutter_fname = outfile_base + '.{0}.stutter_{1}.vcf'.format(region, delta)
            vcf_stutter = get_vcf_writer(stutter_fname)
            mutatant_allele = mutate_str(ref_sequence, repeatunit, delta = delta)
            if delta != 0: # i.e. don't print any lines in the vcf file for the reference allele - it will be a blank vcf.
                ref, alt = trim_indel(ref_sequence, mutatant_allele)
                record = vcf.model._Record(CHROM=chrom, POS=start, ID='.', REF=ref,
                            ALT=[vcf.model._Substitution(alt)],
                            QUAL='.', FILTER='PASS', INFO={'RU':repeatunit},
                            FORMAT='.', sample_indexes=[], samples=None)
                vcf_stutter.write_record(record)
            # write the filename and corresponding stutter probability for use in later pipeline stages
            vcf_probs_writer.write('{0}\t{1}\t{2}\t{3}\n'.format(stutter_fname, prob, bed_out, delta))

if __name__ == '__main__':
    main()
