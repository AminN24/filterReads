#!/usr/bin/env python3

import gzip
import random
import sys
import os
import shutil
import argparse


def read_fastq(r1_stream, r2_stream):
    """ 
        read_fastq() returns two lists. Each list contains the 4 lines of
        information about one of the reads in read pair. The output is in
        binary format.
        --------------
        Line 1: name of the sequence read
        Line 2: sequence
        Line 3: '+'
        Line 4: base quality scores
    """

    r1 = [line for i, line in zip(range(4), r1_stream)]
    r2 = [line for i, line in zip(range(4), r2_stream)]

    return r1, r2


def qual(r1_qual, r2_qual):
    """ 
        qual() returns the quality score for a read pair as the minimum of 
        read 1 and read 2's mean quality scores.

        Inputs:
        -------
        r1_qual: str
            quality scores of read 1 in binary format.

        r2_qual: str
            quality scores of read 2 in binary format.
        -------

        output: float64
            quality score of the read pair.
    """

    r1_q = sum(r1_qual) / len(r1_qual)
    r2_q = sum(r2_qual) / len(r2_qual)
    
    return min(r1_q, r2_q) - 33


def main():

    # ========================================================================
    # ====================== Read inputs from terminal =======================
    # ========================================================================
    parser = argparse.ArgumentParser("Filter read pairs based on quality.")
    parser.add_argument('-1', '--read1', help="read 1 file [required]")
    parser.add_argument('-2', '--read2', help="read 2 file [required]")
    parser.add_argument('-o', '--output',
                    help="output file prefix (path and name) [required]")
    parser.add_argument('-k', type=int, default=5000000,
                    help="number of read pairs to keep [default: %(default)s]")
    parser.add_argument('-v', '--verbose', action='store_true',
                    help="verbose [optional]")
    parser.add_argument('-i', '--interval', type=int, default=200000,
                    help="reporting interval if verbose [default: %(default)s]")

    args = parser.parse_args()

    if args.read1 is None or args.read2 is None or args.output is None:
        parser.print_help()
        sys.exit("Error: missing argument(s)!")

    r1_name = args.read1
    r2_name = args.read2
    out_prefix = args.output
    k = args.k
    rep_interval = args.interval
    verbose = args.verbose
    r1_outpath = f"{out_prefix}_1.subset.fastq.gz"
    r2_outpath = f"{out_prefix}_2.subset.fastq.gz"
    
    # ========================================================================
    # =============== (1) Collect the read pair quality scores ===============
    # ========================================================================
    with gzip.open(r1_name, mode='rb') as r1_stream, \
         gzip.open(r2_name, mode='rb') as r2_stream:

        pair_index = 0
        quals = []
    
        while True:

            # Assign a list to "r1", "r2" containing relevant
            # data for one pair of FASQ reads:
            r1, r2 = read_fastq(r1_stream, r2_stream)

            # Break out of this loop if we've reached the end of file
            if not (r1 and r2):
                print(f"Total number of read pairs: {pair_index}")
                break
        
            pair_index += 1
            quals.append((pair_index, qual(r1[-1].strip(), r2[-1].strip())))

            if verbose and pair_index % rep_interval == 0:
                print(f"Collected quality scores of {pair_index} read pairs.")

    n = len(quals)
    
    # Copying the entire input if subset size is larger than read counts
    if k > n:
        shutil.copyfile(r1_name, r1_outpath)
        shutil.copyfile(r2_name, r2_outpath)
        sys.exit("Subset size is larger than input! Keeping all read pairs.")
    
    # ========================================================================
    # ========== (2) Shuffle for unbiased selection in case of ties ==========
    # ========================================================================
    random.shuffle(quals)

    # ========================================================================
    # ========== (3) Sort read indices based on their quality scores =========
    # ========================================================================
    quals.sort(key=lambda item: item[1])
    
    # Brief report
    if verbose:
        print(f"Original quality score range: {quals[0][1]} - {quals[-1][1]}")
        print(f"Selected quality score range: {quals[(n - k) // 2][1]} - "
              f"{quals[(n + k) // 2 - 1][1]}")

    # ========================================================================
    # ========= (4) Load read indices into a set in order to look up =========
    # ====================== membership in constant time =====================
    # ========================================================================
    index = {pair_index for pair_index, q in quals[(n - k) // 2:(n + k) // 2]}
    
    # ========================================================================
    # ============ (5) Filter input for the selected read pairs ==============
    # ========================================================================
    with gzip.open(r1_name, mode='r') as r1_stream, \
         gzip.open(r2_name, mode='r') as r2_stream, \
         gzip.open(r1_outpath, mode='wb') as r1_record, \
         gzip.open(r2_outpath, mode='wb') as r2_record:

        for pair_index in range(1, n + 1):
        
            r1, r2 = read_fastq(r1_stream, r2_stream)
        
            if pair_index in index:
                r1_record.writelines(r1)
                r2_record.writelines(r2)

            if verbose and pair_index % rep_interval == 0:
                print(f"Filtered {pair_index} read pairs.")


if __name__ == "__main__":
    
    main()
