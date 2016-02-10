#!/usr/bin/env python
""" The ioHCA module regroups the functions linked to file inputs and outputs
"""
from __future__ import print_function
import os, sys, gzip
from Bio import SeqIO

__author__ = "Tristan Bitard-Feildel"
__licence__= "MIT"
__version__ = 0.1
__email__ = "t.bitard.feildel [you know what] uni-muenster.de"
__institute__ = "Institute for Evolution and Biodiversity, Muenster Germany"


def read_multifasta(path, verbose=False):
    """ use Bio.SeqIO to read the fasta file and convert it into a dictionary
    
    Parameter
    ---------
    path : string
        path to the fasta file
        
    Return
    ------
    record_dict : dict
        a dictionary containing the Bio sequence object
    """
    if verbose:
        print("Read fasta inputfile ..")
    if os.path.splitext(path) in [".gz", ".gzip"]:
        with gzip.open(path, 'rt', encoding='utf-8') as handle: #Python3 fix
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        with open(path, "rU") as handle:
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return record_dict
    
    
def write_annotHCA(output, dannotate, sizes, verbose=False):
    """ write annotation output
    
    Parameters
    ----------
    outputf: string
        path to the output file
    dannotate: dict
        the annotation per proteins
    sizes: dict
        protein sizes
    verbose: bool
        print useful stuff
    """
    if verbose:
        print("Writting output of annotation to file {}".format(output))
            
    with open(output, "w") as outf:
        for prot in dannotate:
            outf.write(">{} {}\n".format(prot, sizes[prot]))
            for annotation in dannotate[prot]:
                outf.write("{}\n".format(str(annotation)))

        