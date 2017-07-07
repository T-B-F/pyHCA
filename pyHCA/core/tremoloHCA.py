#!#/usr/bin/env python
""" this is an entirely revised version of tremoloHCA:

NEW features:
    - perform hca segmentation on fasta sequence before if not provided
    - option to select whole sequence or list of domain from hca segmentation
    - iterative search now use hhblits
    - display results by domain arrangement on a text file
    - textfile 2 html is a side utilitary script using tremolo text file result
"""

import os, sys, argparse
from pyHCA.core.annotateHCA import _annotation_aminoacids as segmentation
from pyHCA.core.ioHCA import read_multifasta, write_tremolo_results
from pyHCA.core.external import targets_hhblits, cdd_search, interpro_search
from pyHCA.core.classHCA import Seq

## domains
def read_domainpos(query, positions):
    """ read domain position in tuple format 1,10 20,30
    """
    domains = []
    if positions == None:
        if len(query.seq) == 1:
            # perform segmentation if it's only one sequence
            seg = segmentation(str(query.seq[0]))
            for dom in seg["domain"]:
                domains.append((dom.start, dom.stop))
        else:
            #use orphhca on multiple sequence alignments
            # TODO
            print("Error, hca domain detection not yet implemented", file=sys.stderr)
            sys.exit(1)
    elif positions[0] == "whole":
        # use the whole sequence
        domains = [(0, query.length)]
    else:
        # use user defined positions
        for val in positions:
            start, stop = val.split(",")
            start, stop = int(start)-1, int(stop)
            if stop<=start or start <0 or stop < 1:
                print("Error in start stop values {}".format(val), file=sys.stderr)
            domains.append((start, stop))
    return domains

## targets
def search_domains(query, domains, database, hhblits_evalue, parameters, workdir):
    """ look for targets
    """
    targets = dict()
    alltargetids = set()
    for i, (start, stop) in enumerate(domains):
        #print("domain {}".format(i))
        pathdom = os.path.join(workdir, "dom_{}".format(i))
        if not os.path.isdir(pathdom):
            os.makedirs(pathdom)
        # use sub part of sequence to search for targets
        pathquery = os.path.join(pathdom, "query_{}.fasta".format(i))
        with open(pathquery, "w") as outf:
            for j, name in enumerate(query.name):
                subseq = str(query.seq[j])[start: stop]
                # IMPORTANT: the name of the sequence will be used as input for hhblits
                # a regular expression is set on "Q query_" to catch input name
                outf.write(">query_{} {} {}-{}\n{}\n".format(i, name, start+1, stop, subseq))
        # perform hhblits
        subtargets = targets_hhblits(pathquery, pathdom, database, hhblits_evalue, parameters)
        alltargetids = alltargetids.union(set(subtargets.keys()))
        targets[i] = subtargets
    return targets, list(alltargetids)


## group results by domain arrangements
def group_resda(targets, cddres):
    """ group results per domain arrangements
    """
    groups = dict()
    for querydom in targets:
        groups[querydom] = dict()
        for prot in targets[querydom]:
            if prot in cddres:
                cddres[prot].sort()
                da = ";".join(set([elmnt[2] for elmnt in cddres[prot]]))
                groups[querydom].setdefault(da, list()).append(prot)
            else:
                groups[querydom].setdefault("None", list()).append(prot)
    return groups


def get_cmd():
    """ get command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="inputfasta", 
            help="input fasta file", required=True)
    parser.add_argument("-d", action="store", dest="domains", nargs="+", 
            help="list of domain positions (start and stop inclusive and "
            "separated by comma : -d 1,10 20,30 60,100. If not provided "
            "the search will be performed on each domain found after "
            "segmentation of the input sequence. "
            "To use the whole protein use -d whole.", required=True)
    parser.add_argument("-w", action="store", dest="workdir",
            help="working directory", required=True)
    #parser.add_argument("-a", action="store", dest="annotation", 
            #choices=["CDD", "Interpro"], default="Interpro",
            #help="defined annotation method to use (default=%(default)s)")    
    parser.add_argument("--p2ipr", action="store", dest="p2ipr",
            help="path to the Interpro annotation of UniproKBt proteins, "
                 "gzip format supported.")
                #"If the argument is not specified and '-a Interpro' is set, "
                #"the annotation will be retrieve using web queries of Biomart"
                #" service which will be slower.")
    parser.add_argument("-E", action="store", dest="evalue", 
            help="filter hhblits results by evalue (default=%(default)f)", 
            type=float, default=0.001)
    parser.add_argument("-o", action="store", dest="output", 
            help="output file")
    parser.add_argument("--hhblits-params", action="store", dest="hhblitsparams",
            help="parameters to pass to hhblits, between quotes", default="")
    parser.add_argument("--hhblits-db", action="store", dest="hhblitsdb", 
            help="path to the database to use with hhblits", required=True)
    params = parser.parse_args()
        
    return params

#### MAIN
def main():
    # main tremolo program
    params = get_cmd()

    if not os.path.isdir(params.workdir):
        os.makedirs(params.workdir)

    # read input sequence
    inputquery = read_multifasta(params.inputfasta)
    names, seqs, descrs = list(), list(), list()
    for record in inputquery:
        names.append(inputquery[record].id)
        descrs.append(inputquery[record].description)
        seqs.append(str(inputquery[record].seq))
    query = Seq(names, descrs, seqs, len(seqs[0]))

    # domains? whole sequence? segmentation?
    domains = read_domainpos(query, params.domains)

    # perform search method on each selected parts
    targets, alltargetids = search_domains(query, domains, params.hhblitsdb, params.evalue, params.hhblitsparams, params.workdir)
    if alltargetids == []:
        print("Unable to find any targets with hhblits in database {}".format(params.hhblitsdb), file=sys.stderr)
        print("with parameters {}".format(params.hhblitsparams), file=sys.stderr)
        print("Please try less stringent parameters or a different database", file=sys.stderr)
        with open(params.output, "w") as outf:
            outf.write("# Unable to find any targets with hhblits in database {}\n".format(params.hhblitsdb))
            outf.write("# with parameters {}\n".format(params.hhblitsparams))
            outf.write("# Please try less stringent parameters or a different database\n")

        sys.exit(0)

    #if params.annotation == "CDD":
        # get domain annotation from CDD
        #annotation = cdd_search(alltargetids, params.workdir)
    #else:
        # get domain from Interpro
    annotation = interpro_search(alltargetids, params.workdir, params.p2ipr)

    # group by domain arrangement
    groups = group_resda(targets, annotation)

    # write output
    write_tremolo_results(query, domains, targets, annotation, groups, params.output)

    sys.exit(0)

if __name__ == "__main__":
    main()


