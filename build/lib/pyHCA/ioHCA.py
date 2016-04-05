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

## reorganize results
def flatres(targets, proteins):
    """ flatten results for sorting from a list of proteins
    """
    # flatten dictionary
    #flattargets = list()
    flattargets = dict()
    prot_orders = dict()
    for name in proteins:
        flattargets[name] = dict()
        for hitnum in targets[name]:
            res = targets[name][hitnum]
            if name in prot_orders:
                if res["E-value"] < prot_orders[name]:
                    prot_orders[name] = res["E-value"]
            else:
                prot_orders[name] = res["E-value"]
            flattargets[name][hitnum] = [res["E-value"], res["descr"], 
                res["Probab"], res["Score"], res["Identities"], res["Similarity"], res["Sum_probs"],
                res["Qstart"], res["Qstop"], res["Tstart"], res["Tstop"], 
                res["Qali"], res["Qcons"], res["Tali"], res["Tcons"]]
            #flattargets.append([name, hitnum, res["E-value"], res["descr"], 
                #res["Probab"], res["Score"], res["Identities"], res["Similarity"], res["Sum_probs"],
                #res[d"Qstart"], res["Qstop"], res["Tstart"], res["Tstop"], 
                #res["Qali"], res["Qcons"], res["Tali"], res["Tcons"]])
    flatorders = sorted([(prot_orders[prot], prot) for prot in prot_orders])
    evalues, order = zip(*flatorders) 
    #flattargets.sort()
    return order, flattargets

def orderda(arrangements, targets):
    """ order domain arrangements according to best evalue in the set of hit
    """
    keptda = list()
    for da in arrangements:
        #if da == None:
            #continue
        for prot in arrangements[da]:
            for hitnum in targets[prot]:
                keptda.append((targets[prot][hitnum]["E-value"], da))
    # sort 
    orderedda = list()
    keptda.sort()
    visited = dict()
    for evalue, da in keptda:
        if da not in visited:
            print(evalue, da)
            orderedda.append(da)
            visited[da] = 1
    return orderedda
            
def write_tremolo_results(targets, cddres, groups, output):
    """ write grouped results for domain res
    
    Parameters 
    ---------- 
    targets: dict
        contains for each query domain, the protein and the hits from hhblits
    cddres: dict
        contains the domain annotation
    groups: dict
        group the proteins per domain arrangement
    output: string
        path to the output file

    """
    with open(output, "w") as outf:
        for querydom in groups:
            outf.write("# Domain number {}\n\n".format(querydom))
            # order da depending on best evalue
            orderedda = orderda(groups[querydom], targets[querydom])
            # write all domain arrangement at the beginning an dthe number of proteins
            outf.write("# Domain Domain_arrangement number_of_protein\n")
            for da in orderedda:
                outf.write("INFO\t{}\t{}\t{}\n".format(querydom, da, len(groups[querydom][da])))
            outf.write("\n")
            for da in orderedda:
                # sort prot by evalues
                proteins = groups[querydom][da]
                outf.write("##\n")
                order, flat = flatres(targets[querydom], proteins)
                for prot in order:
                    outf.write(">{}\n".format(prot))
                    if prot in cddres:
                        for start, stop, dom, d_e_val, bitscore, types in cddres[prot]:
                            outf.write("domain\t{}\t{}\t{}\t{}\t{}\t{}\n".format(querydom, dom, start+1, stop, d_e_val, bitscore))
                    else:
                        outf.write("domain\t{}\t{}\n".format(querydom, "None"))
                    for hit in flat[prot]:
                        e_val, descr, prob, score, ident, sim, sprob, qstart, qstop, tstart, tstop, qali, qcons, tali, tcons = flat[prot][hit]
                        outf.write("Hit\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(querydom, hit, e_val, prob, score, ident, sim))
                        outf.write("HitQali\t{}\t{}\t{}\t{}\t{}\n".format(querydom, hit, qstart+1, qstop, qali))
                        outf.write("HitQcon\t{}\t{}\t{}\t{}\t{}\n".format(querydom, hit, qstart+1, qstop, qcons))
                        outf.write("HitTcon\t{}\t{}\t{}\t{}\t{}\n".format(querydom, hit, tstart+1, tstop, tcons))
                        outf.write("HitTali\t{}\t{}\t{}\t{}\t{}\n".format(querydom, hit, tstart+1, tstop, tali))
                        outf.write("//\n")
