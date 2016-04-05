#!/usr/bin/env python
""" this is an entirely revised version of tremoloHCA:

NEW features:
    - perform hca segmentation on fasta sequence before if not provided
    - option to select whole sequence or list of domain from hca segmentation
    - iterative search now use hhblits
    - display results by domain arrangement on a text file
    - textfile 2 html is a side utilitary script using tremolo text file result
"""

import os, sys, argparse
from pyHCA.annotateHCA import _annotation_aminoacids as segmentation
from pyHCA.ioHCA import read_multifasta, write_tremolo_results
from pyHCA.external import targets_hhblits, cdd_search

## domains
def read_domainpos(query, positions):
    """ read domain position in tuple format 1,10 20,30
    """
    domains = []
    if positions == None:
        # perform segmentation
        seg = segmentation(str(query.seq))
        for dom in seg["domain"]:
            domains.append((dom.start, dom.stop))
    elif positions[0] == "whole":
        # use the whole sequence
        domains = [(0, len(query.seq))]
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
def search_domains(query, domains, database, parameters, workdir):
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
        subseq = str(query.seq)[start: stop]
        pathquery = os.path.join(pathdom, "query_{}.fasta".format(i))
        with open(pathquery, "w") as outf:
            outf.write(">query_{} {}-{}\n{}\n".format(i, start+1, stop, subseq))
        # perform hhblits
        subtargets = targets_hhblits(pathquery, pathdom, database, parameters)
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
            

def get_cmd():
    """ get command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", action="store", dest="inputfasta", help="input fasta file", required=True)
    parser.add_argument("-d", action="store", dest="domains", nargs="+", help="list of domain positions (start and stop inclusive and separated by comma : -d 1,10 20,30 60,100. If not provided the search will be performed on each domain found after segmentation of the input sequence. To use the whole protein use -d all.")
    parser.add_argument("-w", action="store", dest="workdir", help="working directory")
    parser.add_argument("-o", action="store", dest="output", help="output file")
    parser.add_argument("--hhblits-params", action="store", dest="hhblitsparams", help="parameters to pass to hhblits, between quotes")
    parser.add_argument("--hhblits-db", action="store", dest="hhblitsdb", help="path to the database to use with hhblits")
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
    if len(inputquery) > 1:
        print("Error, the query should contain only one sequence", file=sys.stderr)
        sys.exit(1)
    for record in inputquery
        query = Seq(record.id, record.descr, str(record.seq))

    # domains? whole sequence? segmentation?
    domains = read_domainpos(query, params.domains)

    # perform search method on each selected parts
    targets, alltargetids = search_domains(query, domains, params.hhblitsdb, params.hhblitsparams, params.workdir)

    # get domain annotation from CDD
    cddres = cdd_search(alltargetids, params.workdir)

    # group by domain arrangement
    groups = group_resda(targets, cddres)

    # write output
    write_tremolo_results(targets, cddres, groups, params.output)

    sys.exit(0)

if __name__ == "__main__":
    main()


