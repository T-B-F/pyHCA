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
            for annotation in dannotate[prot]["domain"]:
                outf.write("{}\n".format(str(annotation)))
            for annotation in dannotate[prot]["cluster"]:
                outf.write("{}\n".format(str(annotation)))

        
def read_annotation(inputfile, formatf):
    """ read domain annotation from pfam or HCA domain file

    Parameters
    ----------
    inputfile: string
        path to the annotation
    formatf: string
        file format either pfam or seghca

    Return
    ------
    annotaton: dict
        the domain annotation
    """
    annotation = dict()
    if formatf == "seghca":
        annotation = read_hcadomain(inputfile)
    elif formatf == "pfam":
        annotation = read_pfamdomain(inputfile)
    else:
        print("Error, no function defined to read domains in format {}".format(formatf), file=sys.stderr)
        sys.exit(1)
    return annotation

def read_hcadomain(inputfile):
    """ read hca domain

    Parameters    
    ---------- 
    inputfile: string
        path to the annotation

    Return
    ------
    annotaton: dict
        the domain annotation
    """
    annotation = dict()
    with open(inputfile) as inf:
        for line in inf:
            if line[0] == ">":
                prot, size = line[1:-1].split()
                annotation[prot] = []
            else:
                tmp = line.split()
                if tmp[0] == "domain":
                    start, stop = int(tmp[1])-1, int(tmp[2])
                    annotation[prot].append((start, stop, tmp[0], "!", None))
    return annotation

def read_pfamdomain(inputfile):
    """ read pfam domain

    Parameters 
    ---------- 
    inputfile: string
        path to the annotation

    Return   
    ------ 
    annotaton: dict 
        the domain annotation
    """  
    annotation = dict()
    with open(inputfile) as inf:
        for line in inf: 
            if line[0] == "#" or line[0] == "\n":
                continue
            tmp = line.split()
            prot = tmp[0]
            name = tmp[1]
            start, stop = int(tmp[2]), int(tmp[3])
            status = "!"
            annotation.setdefault(prot, []).append((start, stop, name, status, None))
    return annotation      
