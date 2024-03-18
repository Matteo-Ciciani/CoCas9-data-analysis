import numpy as np
import pandas as pd
import glob
import fnmatch
import os
import gzip
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
import itertools
import re
from pandas import ExcelWriter
import argparse
    
def generate_raw_PAM_counts(filepaths, targetsites, PAM_length):
    """
    Here, we get all of our relavent PAM sequences from the inputted files
    by searching for the targetsites and looking at the flanking region
    """
    reverse_target_sequences = {targetsite: str(Seq(targetsites[targetsite]).reverse_complement()) for targetsite in targetsites}
    all_pams = {targetsite: [] for targetsite in targetsites}
    # Iterate through each file and collect the PAMs of each sequence
    # Checks both forward and reverse reads
    for filepath in filepaths:
        print('Scanning file: ' + os.path.basename(filepath))
        pams = []
        with gzip.open(filepath, 'rt') as fin:
            for record in SeqIO.parse(fin, 'fastq'):
                seq = str(record.seq)
                for targetsite in targetsites:
                    target_seq = targetsites[targetsite]
                    target = seq.find(targetsites[targetsite])
                    if target > -1:
                        index = target + len(target_seq)
                        all_pams[targetsite].append(seq[index:index + PAM_length])
                    else:
                        target = seq.find(reverse_target_sequences[targetsite])
                        if target > -1:
                            index = target
                            all_pams[targetsite].append(str(Seq(seq[index - PAM_length:index]).reverse_complement()))
    return all_pams
    
def analyze_PAM_depletion_data(filepaths, outdir, targetsites, PAM_length=3):
    """
    Given a directory that contains a given file extension and a target sequence,
    do the entire PAM depletion analysis
    """
    # Make sure that dirnames and target sequences are inputted
    if targetsites is None:
        raise Exception('Please specify a target sequence')
    if PAM_length is None or PAM_length < 3:
        raise Exception('Please enter a valid PAM length')
    if not os.path.exists(outdir):
         os.mkdir(outdir)

    all_pams = generate_raw_PAM_counts(filepaths, targetsites, PAM_length)
    letters = ['A', 'T', 'C', 'G']
    all_counters = {targetsite: Counter(all_pams[targetsite]) for targetsite in targetsites}
    for targetsite in targetsites:
        pams = all_pams[targetsite]
        base_counters = [Counter() for x in range(PAM_length)]
        for pam in pams:
            for i, c in enumerate(pam):
                base_counters[i][c] += 1
        raw_PAM_counts = pd.Series(all_counters[targetsite])
        raw_PAM_counts.sort_values(ascending=False)
        raw_counts_df = pd.DataFrame()
        raw_counts_df['PAM'] = raw_PAM_counts.index
        raw_counts_df['Count'] = raw_PAM_counts.values
        single_base_counts = pd.DataFrame(base_counters)
        single_base_frequencies = single_base_counts.divide(single_base_counts.sum(axis=1).iloc[0])
        # Write counts to excel
        writer = ExcelWriter(os.path.join(outdir, 'counts_' + targetsite + '.xlsx'))
        single_base_counts.to_excel(writer, 'Single Base Counts')
        single_base_frequencies.to_excel(writer, 'Single Base Frequencies')
        raw_counts_df.to_excel(writer, 'Raw PAM Counts')
        writer.save()
        print('Saved excel output for ' + targetsite)

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", type=str, help="one or more fastq files", nargs='+')
    parser.add_argument("-o", '--outdir', type=str, help="The path to the output directory", required=True)
    args = parser.parse_args()
    return args.fastq, args.outdir
        
if __name__ == "__main__":
    file_paths, outdir = arg_parser()
    # Describe the targetsite(s) to search for
    targetsites = {'Site_0': 'GTCGCCCTCGAACTTCACCT'}
    # Run the analysis on the inputted filepaths and targetsite for a given variable nucleotide region length
    analyze_PAM_depletion_data(file_paths, outdir, targetsites, PAM_length=8)
