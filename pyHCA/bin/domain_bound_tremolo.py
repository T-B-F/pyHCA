#!/usr/bin/env python
""" report tremolo results directly matching an annotated protein domain
"""

import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_cmd():
    """ get command line argument
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="tremolores", help="tremolo result file")
    parser.add_argument("-o", action="store", dest="output", help="output plot with covered domains")
    parser.add_argument("--tcov", action="store", dest="min_tcov", type=float, default=0.5, help="minimal coverage of annotated domains to match")
    parser.add_argument("--evalue", action="store", dest="max_evalue", type=float, default=0.001, help="maximal evalue allowed")
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
            if line.startswith("Qdom") and len(tmp) == 4:
                domain = tmp[1]
                start, stop = int(tmp[2])-1, int(tmp[3])
                domains.setdefault(domain, dict())
                domains[domain]["QPos"] = (start, stop)
            elif line.startswith(">"):
                Tname = line[1:].strip()
            elif line.startswith("Qdom") and Tname != None:
                domain = tmp[1]
                Tdomain = tmp[2]
                domains[domain].setdefault(Tname, dict())
                start, stop = int(tmp[3])-1, int(tmp[4])
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
                start, stop = int(tmp[3])-1, int(tmp[4])
                domains[domain][Tname]["Hit"][hitnb]["Qali"] = (start, stop)
            elif tmp[0] == "HitTali" and Tname != None:
                domain = tmp[1]
                hitnb = tmp[2]
                start, stop = int(tmp[3])-1, int(tmp[4])
                domains[domain][Tname]["Hit"][hitnb]["Tali"] = (start, stop)
    return domains

def filter_dommatch(tremolo_res, tcov, tevalue):
    """ find positions matchin a domain
    """
    positions = dict()
    sizes = dict()
    for domain in tremolo_res:
        qdom_start, qdom_stop = tremolo_res[domain]["QPos"]
        sizes[domain] = (qdom_start, qdom_stop)
        del tremolo_res[domain]["QPos"]
        for prot in tremolo_res[domain]:
            target_domains = tremolo_res[domain][prot].get("Tpos", [])
            for hit in tremolo_res[domain][prot]["Hit"]:
                evalue = tremolo_res[domain][prot]["Hit"][hit]["evalue"]
                if evalue < tevalue:
                    thit_start, thit_stop = tremolo_res[domain][prot]["Hit"][hit]["Tali"]
                    qhit_start, qhit_stop = tremolo_res[domain][prot]["Hit"][hit]["Qali"]
                    for tdom_start, tdom_stop, tdom in target_domains:
                        if "/" in tdom and tdom.startswith("IPR"):
                            tdom = tdom.split("/")[0]
                        start = max(thit_start, tdom_start)
                        stop = min(thit_stop, tdom_stop)
                        diff = stop-start
                        if diff > 0:
                            c1 = diff / (tdom_stop-tdom_start)
                            if c1 > tcov:
                                offset = 0
                                if tdom_start > thit_start:
                                    offset = tdom_start - thit_start
                                start = qhit_start + offset
                                
                                offset = thit_stop - thit_start
                                if thit_start < tdom_stop < thit_stop:
                                    offset = tdom_stop - thit_start
                                stop = qhit_start + offset
                                #print(prot, tdom, start, stop, tdom_start, tdom_stop, thit_start, thit_stop)
                                positions.setdefault(domain, dict()).setdefault(tdom, list()).append((start, stop))
                                #print(domain, prot, tdom, tdom_start, tdom_stop, thit_start, thit_stop, qhit_start, qhit_stop, evalue, c1, c2)
    return positions, sizes
    
def main():
    params = get_cmd()
    
    tremolo_res = read_tremolo(params.tremolores)
    
    positions, sizes = filter_dommatch(tremolo_res, params.min_tcov, params.max_evalue)
    
    with PdfPages(params.output) as pdf:
        for query in positions:
            qdom_start, qdom_stop = sizes[query]
            size = qdom_stop - qdom_start
            #print(query, len(positions[query]))
            for domain in positions[query]:
                #print(query, domain)
                coverage = np.zeros(size)
                for start, stop in positions[query][domain]:
                    coverage[start: stop] += 1
                fig, ax = plt.subplots()
                ax.set_title(domain)
                idx = np.arange(size)
                ax.bar(idx+qdom_start, coverage)
                #ax.set_xticks(idx)
                #ax.set_xticklabels((idx+qdom_start)
                pdf.savefig()
                plt.close()
                
                    
    
    sys.exit(0)
    
if __name__ == "__main__":
    main()