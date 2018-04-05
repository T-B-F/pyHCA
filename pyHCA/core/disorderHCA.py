#!/usr/bin/env python
""" compute predicted disorder per residue based on HCA profile
"""

import os, sys, argparse
import numpy as np
from pyHCA import HCA 
from pyHCA.core.seq_util import transform_seq, check_if_msa

def compute_disorder(clusters, seq):
    """ the main function used to compute disorder
    """
    size = len(seq)
    windows = 30

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
    dprofile = dict()
    for i in range(len(seq)):
        sub_Hinclust = Hinclust[i: i+windows].sum()/windows
        sub_Pinclust = Pinclust[i: i+windows].sum()/windows
        sub_outside = outside[i: i+windows].sum()/windows
        logproba, proba = compute_seg_score(sub_outside, sub_Hinclust, sub_Pinclust)
        for j in range(i, i+windows):
            dprofile.setdefault(j, list()).append(logproba)
    lprofile = [sum(dprofile[i])/len(dprofile[i]) for i in range(len(seq))]
    return lprofile

def compute_seg_score(outside, Hinside, Pinside):
    """ compute segment score
    """
    #coef = np.asarray([  0.55207689, -12.35924026,   0.51815519])
    #coef = np.asarray([  0.53857892, -13.48815193,   0.43683726])
    coef = np.asarray([1.69536894, -11.21126274, 1.69447527])
    #intercept = 2.35242882
    #intercept = 0.56405864
    intercept = 1.21383201

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
    parser = argparse.ArgumentParser(prog="{} {}".format(os.path.basename(sys.argv[0]), "draw"))
    parser.add_argument("-i", action="store", dest="fastafile", help="the fasta file", required=True)
    parser.add_argument("-o", action="store", dest="outputfile", help="output file in svg format", required=True)
    parser.add_argument("--verbose", action="store_true", dest="verbose", help="print information")
    params = parser.parse_args()
    
    return params

def main():
    params = get_params()
    
    from pyHCA.core.ioHCA import read_multifasta_it

    with open(params.outputfile, "w") as outf:
        for prot, sequence in read_multifasta_it(params.fastafile):
            sequence = str(sequence.seq)
            seq = transform_seq(sequence)
            hca = HCA(seq=seq)
            clusters = hca.get_clusters()
            disorder = compute_disorder(clusters, seq)
            outf.write(">{}\n".format(prot))
            for i in range(len(seq)):
                outf.write("{} {} {}\n".format(i+1, seq[i], disorder[i]))
                
    sys.exit(0)
    
if __name__ == "__main__":
    main()
