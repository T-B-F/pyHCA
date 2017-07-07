#!/usr/bin/env python

""" The searchHCA module is an ensemble of functions to search a specified 
hydrophobic cluster or a list of  hydrophobic cluster (both borned by start and 
stop positions) on sequence databases.
"""

import os, sys, time, re
import subprocess, shlex
import sqlite3, gzip
debug = True

### CDD

def parse_cddres(rawcdd):
    """ parse raw cdd results
    """
    tmp_cddres = dict()
    cddres = dict()
    for line in rawcdd.split("\n"):
        if line.startswith("Q#"):
            q, types, pssm, begin, end, evalue, score, accession, name, incomplete, superfamily = line.strip().split("\t")

            prot = q.split("-")[-1].strip().split("(")[0][1:] # to manage Warning 
            start, stop = int(begin)-1, int(end)
            evalue = float(evalue)
            bitscore = float(score)
            domname = accession
            if types == "specific" or types == "superfamily":
                tmp_cddres.setdefault(prot, []).append((start, stop, domname.strip(), evalue, bitscore, pssm, types))
    # specific is the best match when overlap, so we keep just specific when there is overlap
    for prot in tmp_cddres:
        cddres[prot] = []
        ldelete = []
        # look for specific and get coords
        for d in tmp_cddres[prot]:
            if d[-1] == 'specific':
                todelete = (d[0], d[1])
                ldelete.append(todelete)
        # if superfamily we search if there is a specific (so the best fit)
        for d in tmp_cddres[prot]:
            if d[-1] == 'superfamily':
                todelete = (d[0], d[1])
                if todelete in ldelete:
                    pass
                    # we only keep if we did not see before in specific
                else:
                    cddres[prot].append(d[:-1])
            # specific 
            else:
                cddres[prot].append(d[:-1])
        cddres[prot].sort()
    return cddres

def cdd_search(listofids, workdir):
    """ get CDD annotation from a list of ids
    """
    cmd = "CDDonline.pl"
    
    pathout = os.path.join(workdir, "list_cdd.out")
    with open(pathout, "w") as outf:
        for ids in listofids:
            outf.write("{}\n".format(ids))
        outf.write("\n")
    # TODO to remove after darkproteome analysis
    #
    # ONLY FOR FASTER DEBUGING
    #
    if os.path.isfile(os.path.join(workdir, "tmp_cdd.out")):
        with open(os.path.join(workdir, "tmp_cdd.out")) as inf:
            raw_cddres = inf.read()
    else:
        try:
            #print("Error unable to find tmp_cdd.out")
            raw_cddres = subprocess.check_output(shlex.split(cmd), input="\n".join(listofids), universal_newlines=True)
        except:
            print("Unable to run CDD search for ids {}".format(" ".join(listofids)), file=sys.stderr)
            sys.exit(1)
        with open(os.path.join(workdir, "tmp_cdd.out"), "w") as outf:
            outf.write(raw_cddres)
    # read cdd res
    cddres = parse_cddres(raw_cddres)
    return cddres

### INTERPRO

def interpro_search(listofids, workdir, path_p2ipr):
    """ get Interpro annotation from a list of ids
    """

    if path_p2ipr:
        return local_interpro_search(listofids, path_p2ipr)
    else:
        #return web_interpro_search(listofids, workdir)
        print("Web search not yet implemented")
        sys.exit(1)

#def read_annotation_sqlite3(uniprotids, path_p2ipr):
    #""" read interpro annotation stored in a sqlite3 database
    #"""

    #annotation = dict()
    #t1 = time.time()
    #conn = apsw.Connection(path_p2ipr)
                
    #cur = conn.cursor()
                        
    #for prot in uniprotids:
        #print(prot)
        #name = uniprotids[prot]
        #cur.execute("select start, stop, interprodom, domain from interpro where protein=?", (prot,))
        #for row in cur:
            #annotation.setdefault(name, list()).append((row[0], row[1], row[2]))

    #cur.close()
    #conn.close()

    #t2 = time.time()
    #print("Done in {}".format(t2-t1))
    #return annotation

def read_annotation(uniprotids, path_p2ipr):
    """ read interpro annotation
    """
    annotation = dict()
    t1 = time.time()
    if path_p2ipr[-2:] == "gz":
        with gzip.open(path_p2ipr, "rt") as inf:
            for line in inf:
                tmp = line.strip().split("\t")
                if tmp[0] in uniprotids:
                    start, stop = int(tmp[4])-1, int(tmp[5])
                    prot = uniprotids[tmp[0]]
                    annotation.setdefault(prot, list()).append((start, stop, tmp[1]+"/"+tmp[3], -1, -1, -1))
                elif len(annotation) == len(uniprotids):
                    return annotation
    else:
        with open(path_p2ipr, "r") as inf:
            for line in inf:
                tmp = line.strip().split("\t")
                if tmp[0] in uniprotids:
                    start, stop = int(tmp[4])-1, int(tmp[5])
                    prot = uniprotids[tmp[0]]
                    annotation.setdefault(prot, list()).append((start, stop, tmp[1]+"/"+tmp[3], -1, -1, -1))
                elif len(annotation) == len(uniprotids):
                    return annotation
    t2 = time.time()
    print("Done in {}".format(t2-t1))
    return annotation

def local_interpro_search(listofids, path_p2ipr):
    """ read local interpro file to look for domain annotations
    """
    annotation = dict()
    uniprotids = dict([(prot.split("|")[1], prot) for prot in listofids])
    cnt = 0
    annotation = read_annotation(uniprotids, path_p2ipr)
    for prot in uniprotids:
        name = uniprotids[prot]
        if name in annotation:
            annotation[name].sort()
        else:
            annotation[name] = list()
    return annotation

### HHBLITS

def filter_hhblits(res, evalue): 
    """ filter targets by only keeping target with at least one hit under <evalue>
    """
    targets = dict()
    for k in res:
        hitnums = list()
        for hitnum in res[k]:
            if res[k][hitnum]["E-value"] < evalue:
                hitnums.append(hitnum)
        if hitnums != list():
            targets[k] = dict()
            for hitnum in hitnums:
                targets[k][hitnum] = dict()
                for param in res[k][hitnum]:
                    targets[k][hitnum][param] = res[k][hitnum][param]
    return targets    

def targets_hhblits(pathquery, workdir, pathdb, cutoff_evalue, parameters):
    """ run and read results of hhblits
    """
    pathout = run_hhblits(pathquery, workdir, pathdb, parameters)
    targets = read_hhblits(pathout)
    new_targets = filter_hhblits(targets, cutoff_evalue)
    return new_targets

def read_hhblits(pathin):
    """ read hhblits results
    """
    alltargets = set()
    targets = dict()
    hitnumber = dict()
    with open(pathin) as inf:
        for line in inf:
            if line[0] == ">":
                tmp = line[1:].split()
                name = tmp[0]
                descr = " ".join(tmp[1:])
                hitnum = 0
                targets.setdefault(name, dict())
                if name in hitnumber:
                    hitnum = hitnumber[name] + 1
                hitnumber[name] = hitnum
                alltargets.add(name)
                targets[name][hitnum] = {"descr": descr, "Tstart":1e10, "Tstop":-1, "Tcons":"", "Tali":"", "Tsize": 0,
                                                 "Qstart":1e10, "Qstop":-1, "Qcons":"", "Qali":"", "Qsize": 0,
                                 "Probab":-1, "E-value":-1, "Score":-1, "Identities":-1, "Similarity":-1, "Sum_probs":-1,
                                }
            elif line.startswith("Probab="):
                #Probab=100.00  E-value=6.2e-40  Score=297.36  Aligned_cols=130  Identities=100%  Similarity=1.267  Sum_probs=129.4
                tmp = line.split()
                for keyval in tmp:
                    key, val = keyval.split("=")
                    if key == "Identities":
                        val = val[:-1]
                    targets[name][hitnum][key] = float(val)
            elif line.startswith("Q ") and line.split()[1] != "Consensus":
                m = re.match("\s*(\d+)\s+([\-\w]*)\s+(\d+)\s\((\d+)\)",line[16:])
                if m:
                    start, ali, stop, size = int(m.group(1))-1, m.group(2), int(m.group(3)), int(m.group(4))
                    targets[name][hitnum]["Qstart"] = min(targets[name][hitnum]["Qstart"], start)
                    targets[name][hitnum]["Qstop"] = max(targets[name][hitnum]["Qstop"], stop)
                    targets[name][hitnum]["Qali"] += ali
                    targets[name][hitnum]["Qsize"] = size
                else:
                    print("FAILED", line, pathin)
            elif line.startswith("Q Consensus"):
                m = re.match("\s*(\d+)\s+([\-\~\w]*)\s+(\d+)",line[16:])
                if m:
                    start, ali, stop = m.group(1), m.group(2), m.group(3)
                    targets[name][hitnum]["Qcons"] += ali
                else:
                    print("FAILED", line, pathin)
            elif line.startswith("T Consensus"):
                m = re.match("\s*(\d+)\s+([\-\~\w]*)\s+(\d+)",line[16:])
                if m:
                    start, ali, stop = m.group(1), m.group(2), m.group(3)
                    targets[name][hitnum]["Tcons"] += ali
                else:
                    print("FAILED", line, pathin)
            elif line.startswith("T "):
                m = re.match("\s*(\d+)\s+([\-\w]*)\s+(\d+)\s\((\d+)\)",line[16:])
                if m:
                    start, ali, stop, size = int(m.group(1)), m.group(2), int(m.group(3)), int(m.group(4))
                    targets[name][hitnum]["Tstart"] = min(targets[name][hitnum]["Tstart"], start)
                    targets[name][hitnum]["Tstop"] = max(targets[name][hitnum]["Tstop"], stop)
                    targets[name][hitnum]["Tali"] += ali
                    targets[name][hitnum]["Tsize"] = size
                elif debug:
                    print("FAILED", line, pathin)
    return targets

def run_hhblits(pathquery, workdir, pathdb, parameters):
    """ run hhblits
    """
    res = dict()
    pathout = os.path.join(workdir, "query_hhblits.hhr")
    #if not  os.path.isfile(pathout):
    pathsco = os.path.join(workdir, "query_hhblits.scores")
    pathlog = os.path.join(workdir, "query_hhblits.log")
    command = "hhblits -i {} -d {} -scores {} -o {} {}".format(pathquery, pathdb, pathsco, pathout, parameters)
    print(command)
    with open(pathlog, "w") as logf:
        try:
            #print("Should no have gone here, hhblits")
            a = subprocess.check_call(shlex.split(command), stdout=logf, stderr=subprocess.STDOUT)
        except:
            print("Unable to run command: {}".format(command), file=sys.stderr)
            sys.exit(1)
    return pathout

