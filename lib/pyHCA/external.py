#!/usr/bin/env python

""" The searchHCA module is an ensemble of functions to search a specified 
hydrophobic cluster or a list of  hydrophobic cluster (both borned by start and 
stop positions) on sequence databases.
"""

import os, sys, subprocess, shlex, re
debug = True

# CDD
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
            domname = name
            if types == "superfamily":
                domname = domname[:-11]
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
    #
    # ONLY FOR FASTER DEBUGING
    #
    #if os.path.isfile(os.path.join(workdir, "tmp_cdd.out")):
        #with open(os.path.join(workdir, "tmp_cdd.out")) as inf:
            #raw_cddres = inf.read()
    #else:
    try:
        raw_cddres = subprocess.check_output(shlex.split(cmd), input="\n".join(listofids), universal_newlines=True)
    except:
        print("Unable to run CDD search for ids {}".format(" ".join(listofids)), file=sys.stderr)
        sys.exit(1)
    #with open(os.path.join(workdir, "tmp_cdd.out"), "w") as outf:
        #outf.write(raw_cddres)
    # read cdd res
    cddres = parse_cddres(raw_cddres)
    return cddres


def targets_hhblits(pathquery, workdir, pathdb, parameters):
    """ run and read results of hhblits
    """
    pathout = run_hhblits(pathquery, workdir, pathdb, parameters)
    targets = read_hhblits(pathout)
    return targets

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
                targets[name][hitnum] = {"descr": descr, "Tstart":1e10, "Tstop":-1, "Tcons":"", "Tali":"",
                                                 "Qstart":1e10, "Qstop":-1, "Qcons":"", "Qali":"",
                                 "Probab":-1, "E-value":-1, "Score":-1, "Identities":-1, "Similarity":-1, "Sum_probs":-1,
                                }
            elif line.startswith("Probab="):
                tmp = line.split()
                for keyval in tmp:
                    key, val = keyval.split("=")
                    if key == "Identities":
                        val = val[:-1]
                    targets[name][hitnum][key] = float(val)
            elif line.startswith("Q query_"):
                m = re.match("\s*(\d+)\s+([\-\w]*)\s+(\d+)",line[16:])
                if m:
                    start, ali, stop = int(m.group(1))-1, m.group(2), int(m.group(3))
                    targets[name][hitnum]["Qstart"] = min(targets[name][hitnum]["Qstart"], start)
                    targets[name][hitnum]["Qstop"] = max(targets[name][hitnum]["Qstop"], stop)
                    targets[name][hitnum]["Qali"] += ali
                elif debug:
                    print("FAILED", line)
            elif line.startswith("Q Consensus"):
                m = re.match("\s*(\d+)\s+([\-\~\w]*)\s+(\d+)",line[16:])
                if m:
                    start, ali, stop = m.group(1), m.group(2), m.group(3)
                    targets[name][hitnum]["Qcons"] += ali
                elif debug:
                    print("FAILED", line)
            elif line.startswith("T Consensus"):
                m = re.match("\s*(\d+)\s+([\-\~\w]*)\s+(\d+)",line[16:])
                if m:
                    start, ali, stop = m.group(1), m.group(2), m.group(3)
                    targets[name][hitnum]["Tcons"] += ali
                elif debug:
                    print("FAILED", line)
            elif line.startswith("T "):
                m = re.match("\s*(\d+)\s+([\-\w]*)\s+(\d+)",line[16:])
                if m:
                    start, ali, stop = int(m.group(1)), m.group(2), int(m.group(3))
                    targets[name][hitnum]["Tstart"] = min(targets[name][hitnum]["Tstart"], start)
                    targets[name][hitnum]["Tstop"] = max(targets[name][hitnum]["Tstop"], stop)
                    targets[name][hitnum]["Tali"] += ali
                elif debug:
                    print("FAILED", line)
    return targets

def run_hhblits(pathquery, workdir, pathdb, parameters):
    """ run hhblits
    """
    res = dict()
    pathout = os.path.join(workdir, "query_hhblits.hhr")
    pathsco = os.path.join(workdir, "query_hhblits.scores")
    pathlog = os.path.join(workdir, "query_hhblits.log")
    command = "hhblits -i {} -d {} -scores {} -o {} {}".format(pathquery, pathdb, pathsco, pathout, parameters)
    with open(pathlog, "w") as logf:
       try:
           a = subprocess.check_call(shlex.split(command), stdout=logf, stderr=subprocess.STDOUT)
       except:
           print("Unable to run command:", file=sys.stderr)
           print(command, file=sys.stderr)
           sys.exit(1)
    return pathout
