#!/usr/bin/env python
""" make spaced seed based on HCA
"""

#from pyHCA.core.ioHCA import read_singlefasta
from pyHCA.core.ioHCA import read_multifasta_it

from pyHCA.core.HCA import HCA

import os, sys, argparse

def read_singlefasta(path, verbose=False):
    """ use Bio.SeqIO to read the fasta file and assert that it only contains 
    one protein
    """
    cnt = 0
    query, sequence = None, None
    for record, seq in read_multifasta_it(path, verbose=False):
        if cnt > 0:
            raise RuntimeError("Error, multiple sequences found in file {}".format(path))        
        query = record
        sequence = str(seq.seq)
        cnt += 1
    return query, sequence

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="inputfasta")
    parser.add_argument("-o", action="store", dest="outputfile")
    params = parser.parse_args()
    return params


def main():
    params = get_cmd()
    
    query, input_seq = read_singlefasta(params.inputfasta)
    hca = HCA(seq=input_seq.upper(), querynames=query)
    clusters = hca.get_clusters(prot=query)
    
    with open(params.outputfile, "w") as outf:
        for clust in clusters:
            if len(clust.hydro_cluster) > 2:
                outf.write("{}\t{}\n".format(query, "\t".join(str(clust).split()[1:])))

    sys.exit(0)

if __name__ == "__main__":
    main()