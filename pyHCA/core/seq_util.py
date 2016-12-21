#!/usr/bin/env python
""" sequences handling functions
"""

import os, sys, argparse, string
from Bio import Seq

__author__ = "Tristan Bitard-Feildel, Guillem Faure"
__licence__= "MIT"
__version__ = 0.1
__email__ = "t.bitard.feildel [you know what] uni-muenster.de"
__institute__ = "Institute for Evolution and Biodiversity, Muenster Germany"

# problem of the SeqIO module: not memory efficient

__all__ = ["itercodon", "six_frames", "transform_seq"]

def transform_seq(seq):
    return seq.replace("*", "").replace("-", "").replace("?", "").replace("!","")

def itercodon(seq, frame, offset, table, reverse=False):
    stop = 0
    if not reverse:
        for i in xrange(frame, len(seq)-offset, 3):
            subseq = str(seq.seq)[i:i+3]
            assert(len(subseq)%3==0),(str(seq))
            aa = Seq.translate(subseq, table)
            yield i, aa
        if i+3 != len(seq):
            subseq = seq[i+3:] + "N"*(3-offset)
            assert(len(subseq)%3==0)
            aa = Seq.translate(subseq, table)
            yield i, aa
    else:
        for i in xrange(len(seq), offset, -3):
            # the reverse complement
            subseq = Seq.reverse_complement(str(seq.seq)[i-3:i])
            assert(len(subseq)%3==0)
            aa = Seq.translate(subseq, table)
            yield i, aa
        if offset:
            subseq = Seq.reverse_complement("N"*(3-offset) + str(seq.seq)[:offset])
            assert(len(subseq)%3==0)
            aa = Seq.translate(subseq, table)
            yield i, aa
    
#modulo = 0
#frame 0 0
#frame 1 -2
#frame 2 -1

#modulo = 1
#frame 0 -1
#frame 1  0
#frame 2 -2 

#modulo = 2
#frame 0 -2
#frame 1 -1
#frame 2  0

def six_frames(seq, table=1):
    """ Return the six frame translation of Seq protein object
    """
    offset = len(seq) % 3
    if len(seq) % 3 == 0:
        offset = [0, 2, 1]
    elif len(seq) % 3 == 1:
        offset = [1, 0, 2]
    else:
        offset = [2, 1, 0]
    start = 0
    reverse, prev = False, False
    for strand in [1, -1]:
        if strand < 0:
            reverse = True
        for frame in range(3):
            subprot = ""
            for nuc_idx, aa in itercodon(seq, frame, offset[frame], table, reverse):
                #if frame == 1 and strand == 1:
                    #print (start, nuc_idx, aa)
                if aa == "*" and subprot:
                    yield strand, frame, start, subprot
                    subprot = ""
                    prev = True
                elif prev:
                    start = nuc_idx
                    prev = False
                else:
                    subprot += aa
            if subprot:
                yield strand, frame, start, subprot
            