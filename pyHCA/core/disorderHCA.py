#!/usr/bin/env python
""" compute predicted disorder per residue based on HCA profile
"""

import os, sys, argparse
import numpy as np
from pyHCA import HCA 
from pyHCA.core.seq_util import transform_seq, check_if_msa

def compute_disorder(clusters, seq, windows=30, 
                     norm=False, intercept=None, coef=None):
    """ the main function used to compute disorder
    """
    size = len(seq)

    Hinclust = np.zeros(size) 
    Pinclust = np.zeros(size) 
    outside = np.zeros(size) 
    outside.fill(1)  
    for clust in clusters:
        if len(clust.hydro_cluster) > 2:
            outside[clust.start: clust.stop] = 0 
            for i in range(len(clust.hydro_cluster)):
                if clust.hydro_cluster[i] == 1:
                    Hinclust[clust.start+i] = 1
                else:
                    Pinclust[clust.start+i] = 1
    if len(seq) > windows:
        dprofile = dict()
        for i in range(len(seq)-windows+1):
            sub_Hinclust = Hinclust[i: i+windows].sum()
            sub_Pinclust = Pinclust[i: i+windows].sum()
            sub_outside = outside[i: i+windows].sum()
            if norm:
                sub_Hinclust /= windows
                sub_Pinclust /= windows
                sub_outside  /= windows
            logproba, proba = compute_seg_score(sub_outside, sub_Hinclust, sub_Pinclust, 
                                                windows, intercept, coef)
            for j in range(i, i+windows):
                dprofile.setdefault(j, list()).append(logproba)
    else:
        return dict()
        #dprofile = dict()
        #for i in range(-windows+1, len(seq)):
        #    sub_size = i+min(windows, len(seq))
        #    sub_Hinclust = Hinclust[max(0, i): i+min(len(seq), windows)].sum()/sub_size
        #    sub_Pinclust = Pinclust[max(0, i): i+min(len(seq), windows)].sum()/sub_size
        #    sub_outside = outside[max(0, i): i+min(len(seq), windows)].sum()/sub_size
        #    logproba, proba = compute_seg_score(sub_outside, sub_Hinclust, sub_Pinclust, windows)
        #    for j in range(i, i+windows):
        #        dprofile.setdefault(j, list()).append(logproba)
    lprofile = [sum(dprofile[i])/len(dprofile[i]) for i in range(len(seq))]
    return lprofile

def compute_seg_score(outside, Hinside, Pinside, windows, intercept=None, coef=None):
    """ compute segment score
    """
    parameters = {
        #windows: (intercept, coef)
        10: (0.82510422, np.asarray([ 0.04843165, -0.43583164,  0.        ])),
        15: (1.3530656, np.asarray([ 0.43132414, -6.82680033,  0.33753246])),
        20: (-2.49828023, np.asarray([ 0.24618441, -0.23579668,  0.2280128 ])),
        25: (0.42444814, np.asarray([ 0.08732703, -0.36992488,  0.0930069 ])),
        #30: (1.19777441, np.asarray([  1.66817336, -11.11173207,   1.60691731])),
        30: (-12.36880161, np.asarray([ 24.95909262, -36.28374948,  15.06856888])),
        35: (-2.99280193, np.asarray([ 0.18838551, -0.23499759, 0.19102656])),

    }
    if intercept is None:
        intercept, coef = parameters[windows] 
    x = np.asarray([outside, Hinside, Pinside], dtype=float)
    V = x * coef
    U = V.sum() + intercept
    U *= -1
    A = np.exp(U)+1
    logP = -np.log(A)
    P = np.reciprocal(A)
    return logP, P


def get_params():
    """ get command line ArgumentParser
    """
    parser = argparse.ArgumentParser(prog="{} {}".format(os.path.basename(sys.argv[0]), "disorder"))
    parser.add_argument("-i", action="store", dest="fastafile", help="the fasta file", required=True)
    parser.add_argument("-w", action="store", dest="windows", help="the windows size", 
                        choices=[10, 15, 20, 25, 30, 35], type=int, required=True)
    parser.add_argument("-o", action="store", dest="outputfile", help="output file in svg format", required=True)
    parser.add_argument("--normed", action="store_true", dest="norm", help="norm data")
    parser.add_argument("--intercept", action="store", dest="intercept", help="intercept")
    parser.add_argument("--coefs", action="store", dest="coef", help="coefficients", nargs="+")
    parser.add_argument("--verbose", action="store_true", dest="verbose", help="print information")
    params = parser.parse_args()
    
    return params

def main():
    params = get_params()
    
    from pyHCA.core.ioHCA import read_multifasta_it

    if params.intercept is not None:
        params.intercept = float(params.intercept)
        params.coef = np.asarray([float(v) for v in params.coef])

    with open(params.outputfile, "w") as outf:
        for prot, sequence in read_multifasta_it(params.fastafile):
            sequence = str(sequence.seq)
            seq = transform_seq(sequence)
            hca = HCA(seq=seq)
            clusters = hca.get_clusters()
            disorder = compute_disorder(clusters, seq, params.windows, params.norm, params.intercept, params.coef)
            if disorder:
                outf.write(">{}\n".format(prot))
                for i in range(len(seq)):
                    outf.write("{} {} {}\n".format(i+1, seq[i], disorder[i]))
            else:
                outf.write(">{} Error, sequence too small ({}, expected {})\n".format(prot, len(seq), params.windows))
                
    sys.exit(0)
    
if __name__ == "__main__":
    main()

