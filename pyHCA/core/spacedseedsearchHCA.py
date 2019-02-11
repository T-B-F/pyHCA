#!/usr/bin/env python
""" spaced seed search based on HCA

Seeds:
10001
110101101
1001001
...

Sequences
00001001010010000101010010101010101000101001000001 ...
01101001010100101001001010110100000110000001001000 ...

Keep sequence_i if:
- number of seeds found (conserved order) above threshold
- percentage of seeds found (/ sum length of seeds) above threshold


Options:
- allow partial seed match
 0000100010100001 ith seed 110101
      110101
 -------***------
- allow seed 0 to match a 1 in sequence
- change thresholds

"""

import os, sys, argparse
from pyHCA.core.HCA import HCA

def get_cmd():
    params = argparse.ArgumentParser()
    params.add_argument("-i", action="store", dest="inputfasta")
    params.add_argument("-o", action="store", dest="outputfile")
    params.add_argument("-w", action="store", dest="workdir")
    params.add_argument("-s", action="store", dest="seedfile")
    params.add_argument("-d", action="store", dest="seqdatabase")
    params.add_argument("--options", action="store", dest="options", 
        help="options to pass to jackhmmer")
    parser = params.parse_args()
    return parser

def read_seed(path):
    """ read spaced seeds / hydrophobic clusters and store their order 
    and positions
    """
    # TODO look for an efficient data structure
    seeds = dict()
    with open(path) as inf:
        for line in inf:
            tmp = line.split()
            start = int(tmp[1])-1
            seed = tmp[3]
            seeds.setdefault(seed, list()).append(start)
    return seeds

def spaced_seed_seq_similarity(sequence, spaced_seed):
    """ search spaced seed in a sequence and return the number of spaced 
    seed found along with coverage (to weight short seeds)
    """
    # TODO look for an efficient data structure and algorithm
    # this part should be fast
    hca = HCA(seq=sequence, querynames="query")
    seqbin = hca.get_seqbin()
    starts = [i for i in range(len(seqbin)) if seqbin[i] == 1]
    
    

def filter_sequences(sequences, spaced_seeds, found_threshold=0.5, cov_threshold=0.5,)
    """ filter a list of sequences based on the maximum number of spaced seed found
    in the correct order and their coverage
    """
    kept = list()
    for i, seq in enumerate(sequences):
        perc_found, cov = spaced_seed_seq_similarity(seq, spaced_seeds)
        if perc_found > found_threshold and cov > cov_threshold:
            kept.append(i)
    return kept

def main():
    params = get_cmd()
    
    seeds = read_seed(params.seedfile)
    
    # prepare and run jackhmmer
    
        sequences = iteration()
        kept = filter_sequences(sequences, seeds)
        # write kept sequences for next iterations
    

    sys.exit(0)

if __name__ == "__main__":
    main()