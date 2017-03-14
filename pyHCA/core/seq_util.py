#!/usr/bin/env python
""" sequences handling functions
"""

import os, sys, argparse, string
from Bio import Seq
from Bio.SubsMat import MatrixInfo

__author__ = "Tristan Bitard-Feildel, Guillem Faure"
__licence__= "MIT"
__version__ = 0.1
__email__ = "t.bitard.feildel [you know what] uni-muenster.de"
__institute__ = "Institute for Evolution and Biodiversity, Muenster Germany"

# problem of the SeqIO module: not memory efficient

__all__ = ["itercodon", "six_frames", "transform_seq"]

def transform_seq(seq):
    return seq.replace("*", "").replace("-", "").replace("?", "").replace("!","").replace(".", "")

def check_seq_char(seq):
    illegal_char = set(["*", "-", "?", "!", "."])
    for c in seq:
        if c in illegal_char:
            return True
    return False

def is_msa(sequences):
    if isinstance(sequences, dict):
        for rec in sequences:
            seq = sequences[rec]
            if check_seq_char(seq):
                return True
    elif isinstance(sequences, list):
        for seq in sequences:
            if check_seq_char(seq):
                return True
    elif isinstance(sequences, str):
        seq = sequences[:]
        if check_seq_char(seq):
            return True
    else:
        raise ValueError("Unknown argument type passed to is_msa(), {}").format(type(sequences))
    return False
    

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
    
def compute_conserved_positions(dfasta, dmsa, score_type=0, matrix_name="blosum62"):
    """ compute conservation of a column relative to a sequence position of a msa
    score_type, 0: identity score, 1: binarized similarity (1 if sim > 0 else 0)
    """ 
    
    matrix = getattr(MatrixInfo, matrix_name)
    items = list(matrix.items())
    matrix.update(((b,a),val) for (a,b),val in items)
    
    dconserv_per_prot = dict()
    records = list()
    seq_idx = dict()
    for rec in dfasta:
        records.append(rec)
        seq_idx[rec] = 0
        dconserv_per_prot[rec] = [0] * len(dfasta[rec])
        
    nb_seq =  len(records)
    if nb_seq > 0:
        nb_cols = len(dmsa[records[0]])
        for c in range(nb_cols):
            #dconserv[c] = 0.0
            for k in range(len(records)-1):
                i = seq_idx[records[k]]
                seqk = dmsa[records[k]]
                if seqk[c] != "-":
                    # no gap in sequence k
                    for l in range(k+1, len(records)):
                        j = seq_idx[records[l]]
                        seql = dmsa[records[l]]
                        if seql[c] != "-":
                            if score_type == 0:
                                score = 1 if (seqk[c] == seql[c]) else 0
                            else: # score_type == 1:
                                score = 1 if matrix[(seqk[c], seql[c])] > 0 else 0
                            dconserv_per_prot[records[k]][i] += score
                            dconserv_per_prot[records[l]][j] += score
            # update indexe position of each sequence
            for rec in records:
                if dmsa[rec][c] != "-":
                    seq_idx[rec] += 1
        # normalize score by number of sequence
        for rec in records:
            for i in range(len(dconserv_per_prot[rec])):
                dconserv_per_prot[rec][i] /= (nb_seq-1)
    return dconserv_per_prot
