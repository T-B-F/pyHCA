#!/usr/bin/env python
""" compute predicted disorder per residue based on HCA profile
"""

import os, sys, argparse
import numpy as np
try:
    from sklearn.externals import joblib
except:
    import joblib
from pyHCA import HCA 
from pyHCA.core.seq_util import transform_seq, check_if_msa

def prepare_sequence(seq):
    """ compute score for a given sequence
    """
    size = len(seq)
    hcaprot = HCA(seq=seq)
    hcaclusters = hcaprot.get_clusters()
    hcadomains  = hcaprot.get_domains()
    return  hcadomains, hcaclusters

polar = ["R", "K", "N", "D", "E", "Q"]
def compute_features2(seq, domains, clusters, AA_sorted):
    """ list of features:
        - amino acid type
        - domain size
        - cluster size
        - minimal distance to an other cluster (or termini)
        - cluster density
        - number of hydrophobic in -8 i
        - number of hydrophobic in +8 i
        - number of polar in -8 i
        - number of polar in +8 i
        - number of polar in cluster
        - number of hydrophobic in cluster
        - cluster propensity to H,E, h,e, d/c?
    n = 10
    """
    n = 9
    X = np.zeros((len(seq), n))
    domain_size = dict()
    for dom_id, dom in enumerate(domains):
        for i in range(dom.start, dom.stop):
            domain_size[i] = (dom_id, dom.stop - dom.start + 1)

    cluster_size = dict()
    cluster_hydro = dict()
    cluster_pos = np.zeros(len(seq))
    number_hydro_clust = dict()
    cluster_positions = list()
    for cls_id, cls in enumerate(clusters):
        cluster_positions.append((cls.start, cls.stop))
        number_hydro_clust[cls_id] = 0
        k = 0
        for i in range(cls.start, cls.stop):
            cluster_pos[cls.start: cls.stop] = 1
            cluster_size[i] = (cls_id, cls.stop - cls.start + 1)
            if cls.hydro_cluster[k] == 1:
                number_hydro_clust[cls_id] += 1
            k+=1
    polar_seq = np.zeros(len(seq))
    for i in range(len(seq)):
        if seq[i] in polar:
            polar_seq[i] = 1
            
    cluster_positions.sort()
    for i in range(len(seq)):

        dom_id, dom_size = domain_size.get(i, (-1, 0))
        cls_id, cls_size = cluster_size.get(i, (-1, 0))             
        number_hydro = 0

        if cls_id >= 0:
            number_hydro = number_hydro_clust[cls_id] / cls_size
            cluster = cluster_positions[cls_id]
            dist_to_beg = cluster[0]
            dist_to_end = len(seq)-cluster[1] +1
            min_cluster_distance = None
            if cls_id > 0:
                prev_cluster = cluster_positions[cls_id-1]
                d_prev = cluster[0] - prev_cluster[1]
                dist_to_beg = min(d_prev, dist_to_beg)
            if cls_id < len(cluster_positions) - 1:
                next_cluster = cluster_positions[cls_id+1]
                d_next = next_cluster[1] - cluster[1]
                dist_to_end = min(d_next, dist_to_end)
            min_cluster_distance = min(dist_to_beg, dist_to_end)
        else:
            dist_to_beg = i
            dist_to_end = len(seq) - i + 1
            min_cluster_distance = min(dist_to_beg, dist_to_end)
            for j in range(len(cluster_positions)):
                d1 = np.abs(cluster_positions[j][0] - i)
                d2 = np.abs(cluster_positions[j][1] - i)
                d = min(d1, d2)
                min_cluster_distance = min(min_cluster_distance, d)
                
        hydro_res_in_cls = sum(cluster_pos[max(0, i-4): min(i+5, len(seq)+1)])
        hdensity1 = hydro_res_in_cls / (min(i+5, len(seq)+1)-max(0, i-4))
        hydro_res_in_cls = sum(cluster_pos[max(0, i-8): min(i+9, len(seq)+1)])
        hdensity2 = hydro_res_in_cls / (min(i+9, len(seq)+1)-max(0, i-8))

        polar_res_in_cls = sum(polar_seq[max(0, i-4): min(i+5, len(seq)+1)])
        pdensity1 = polar_res_in_cls / (min(i+5, len(seq)+1)-max(0, i-4))
        polar_res_in_cls = sum(polar_seq[max(0, i-8): min(i+9, len(seq)+1)])
        pdensity2 = polar_res_in_cls / (min(i+9, len(seq)+1)-max(0, i-8))

        aa = seq[i].upper()
        if aa not in AA_sorted:
            aa = "X"
        aa_int = AA_sorted[aa]
        X[i][0] = aa_int
        X[i][1] = dom_size
        X[i][2] = cls_size
        X[i][3] = min_cluster_distance
        X[i][4] = hdensity1
        X[i][5] = hdensity2
        X[i][6] = number_hydro
        X[i][7] = pdensity1
        X[i][8] = pdensity2
    return X

def compute_features3(seq, domains, clusters, AA_sorted):
    features = compute_features2(seq, domains, clusters, AA_sorted)
    features = features[:, [0, 1, 4, 5, 6, 8]]
    return features

def get_params():
    """ get command line ArgumentParser
    """
    parser = argparse.ArgumentParser(prog="{} {}".format(os.path.basename(sys.argv[0]), "disorder"))
    parser.add_argument("-i", action="store", dest="fastafile", help="the fasta file", required=True)
    parser.add_argument("-o", action="store", dest="outputfile", help="output file in svg format", required=True)
    parser.add_argument("-m", action="store", dest="model", help="model to use", required=True)
    parser.add_argument("--verbose", action="store_true", dest="verbose", help="print information")
    params = parser.parse_args()
    
    return params

def main():
    params = get_params()
    
    from pyHCA.core.ioHCA import read_multifasta_it

    AA_sorted = dict()
    AA1 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] 
    for i, aa in enumerate(AA1):
        AA_sorted[aa] = i+1
    if "X" not in AA_sorted:
        AA_sorted["X"] = 0
                
    compute_features = compute_features3 # all_compute_features[params.method]

    scaler_model = joblib.load(params.model)
    rbs_scaler = scaler_model["scaler"]
    trained_clf = scaler_model["model"]

    with open(params.outputfile, "w") as outf:
        for prot, sequence in read_multifasta_it(params.fastafile):
            seq = str(sequence.seq).upper()
            domains, clusters = prepare_sequence(seq)
            features = compute_features(seq, domains, clusters, AA_sorted)
            features = rbs_scaler.transform(features)
            probas = trained_clf.predict_proba(features)

            outf.write(">{}\n".format(prot))
            for i in range(len(probas)):
                outf.write("{} {} {}\n".format(i+1, seq[i],  probas[i, 0]))
    sys.exit(0)
    
if __name__ == "__main__":
    main()

