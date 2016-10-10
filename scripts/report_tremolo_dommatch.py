#!/usr/bin/env python
""" report tremolo results directly matching an annotated protein domain
"""

import os, sys, argparse

def get_cmd():
    """ get command line argument
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="tremolores")
    parser.add_argument("--qcov", action="store", dest="min_qcov", type=float, default=0.8)
    parser.add_argument("--tcov", action="store", dest="min_tcov", type=float, default=0.5)
    params = parser.parse_args()
    return params

def read_tremolo(path):
    """ read tremolo domain results
    """
    domains = dict()
    Tname = None
    with open(path) as inf :
        for line in inf:
            #print(line)
            if line[0] == "\n" or line[0] == "#":
                Tname = None
                continue
            tmp = line.strip().split("\t")
            if line.startswith("Qdom") and len(tmp) == 3:
                domain = tmp[0].split()[1]
                start, stop = int(tmp[1]), int(tmp[2])
                domains.setdefault(domain, dict())
                domains[domain]["QPos"] = (start, stop)
            elif line.startswith(">"):
                Tname = line[1:].strip()
            elif line.startswith("Qdom") and Tname != None:
                domain = tmp[0].split()[1]
                Tdomain = tmp[1]
                domains[domain].setdefault(Tname, dict())
                start, stop = int(tmp[2]), int(tmp[3])
                domains[domain][Tname].setdefault("Tpos", list()).append((start, stop, Tdomain))
            elif tmp[0] == "Hit" and Tname != None:
                domain = tmp[1]
                domains[domain].setdefault(Tname, dict())
                hitnb = tmp[2]
                evalue = float(tmp[3])
                domains[domain][Tname].setdefault("Hit", dict()).setdefault(hitnb, dict())
                domains[domain][Tname]["Hit"][hitnb]["evalue"] = evalue
            elif tmp[0] == "HitQali" and Tname != None:
                domain = tmp[1]
                hitnb = tmp[2]
                start, stop = int(tmp[3]), int(tmp[4])
                domains[domain][Tname]["Hit"][hitnb]["Qali"] = (start, stop)
            elif tmp[0] == "HitTali" and Tname != None:
                domain = tmp[1]
                hitnb = tmp[2]
                start, stop = int(tmp[3]), int(tmp[4])
                domains[domain][Tname]["Hit"][hitnb]["Tali"] = (start, stop)
    return domains

def find_dommatch(tremolo_res, qcov, tcov):
    """ find positions matchin a domain
    """
    for domain in tremolo_res:
        qdom_start, qdom_stop = tremolo_res[domain]["QPos"]
        del tremolo_res[domain]["QPos"]
        for prot in tremolo_res[domain]:
            target_domains = tremolo_res[domain][prot].get("Tpos", [])
            for hit in tremolo_res[domain][prot]["Hit"]:
                evalue = tremolo_res[domain][prot]["Hit"][hit]["evalue"]
                if evalue < 0.001:
                    thit_start, thit_stop = tremolo_res[domain][prot]["Hit"][hit]["Tali"]
                    qhit_start, qhit_stop = tremolo_res[domain][prot]["Hit"][hit]["Qali"]
                    for tdom_start, tdom_stop, tdom in target_domains:
                        start = max(thit_start, tdom_start)
                        stop = min(thit_stop, tdom_stop)
                        diff = stop-start
                        if diff > 0:
                            c1 = diff / (tdom_stop-tdom_start)
                            c2 = diff / (qdom_stop-qdom_start)
                            if c1 > tcov and c2 > qcov:
                                print(domain, prot, tdom, tdom_start, tdom_stop, thit_start, thit_stop, qhit_start, qhit_stop, evalue, c1, c2)
    
def main():
    params = get_cmd()
    
    tremolo_res = read_tremolo(params.tremolores)
    
    find_match = find_dommatch(tremolo_res, params.min_qcov, params.min_tcov)
    
    sys.exit(0)
    
if __name__ == "__main__":
    main()