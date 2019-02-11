#!/usr/bin/env python
""" make spaced seed based on HCA
"""

#from pyHCA.core.ioHCA import read_singlefasta
from pyHCA.core.ioHCA import read_multifasta
from pyHCA.core.HCA import HCA
import numpy as np
import os, sys, argparse

def check_msa(dfasta):
    """ check that sequences of msa are all the same length
    """
    size = None
    for q in dfasta:
        if size is None:
            size = len(str(dfasta[q].seq))
        elif size != len(str(dfasta[q].seq)):
            raise ValueError("Different sequence lengths detected in MSA")

def insert_gaps(gapped_seq, clusters):
    """ insert gaps inside hydrophobic clusters (computed from ungapped
    sequences of an MSA

    Arguments:
    ----------
    - gapped_seq: Seq
        a Biopython sequence from a MSA
    - clusters:
        hydrophobic clusters computed from the corresponding ungapped sequence

    Returns:
    --------
    - new_clusters: list
        a list of (start, stop, hydrophobic clusters) of which indexes have
        been offseted to take into account gap positions and -1 have been
        added to hydrophobic clusters
    """
    offset = 0
    index = 0
    ungapped_seq = gapped_seq.replace("-", "")
    new_clusters = list()
    for cl in clusters:
        start = cl.start
        stop = cl.stop
        hclst = cl.hydro_cluster
        while index < start:
            if gapped_seq[index+offset] == "-":
                offset += 1
            else:
                index += 1
        new_start = cl.start + offset    
        new_hclst = []
        n = 0
        while index+n < stop:
            if gapped_seq[index+offset+n] == "-":
                offset += 1
                new_hclst.append(-1)
            else:
                new_hclst.append(hclst[n])
                n += 1   
        index += n
        new_stop = cl.stop + offset
        new_clusters.append((new_start, new_stop, new_hclst))
    return new_clusters

def map_clusters_ali(dfasta, clusters):
    """ map clusters of unaligned sequences to MSA
    
    Arguments:
    ----------
    - dfasta: dict
        dictionary of aligned protein sequences
    - clusters: dict
        dictionary of HCA clusters (unaligned sequences)

    Returns:
    --------
    - clusters_pos: dict
        for each protein an array of MSA length with value 
        corresponding to cluster index label
    - clusters_ali: dict
        for each protein a list of clusters (start, stop, array of 
        hydro clusters) with gap inserted in original hydrophobic
        clusters (gaps are labelled as -1)
    """
    clusters_ali = dict()
    clusters_pos = dict()
    for prot in dfasta:
        clusters_with_gaps = insert_gaps(str(dfasta[prot].seq),
                                         clusters[prot])
        ar_0 = np.zeros(len(str(dfasta[prot].seq)))
        for i, cl in enumerate(clusters_with_gaps):
            ar_0[cl[0]: cl[1]] = i + 1
        clusters_pos[prot] = ar_0
        clusters_ali[prot] = clusters_with_gaps
    return clusters_pos, clusters_ali

def cluster_conservation(clusters, size):
    """ compute conservation of hydrophobe/hydrophile/gap characters
    between hydrophobic clusters of MSA

    Arguments:
    ----------
    - clusters: dict
        keys are proteins, values are list of (start, stop, hydrophobic
        clsters) with gaps and offseted
    - size: int
        the msa size
    
    Returns:
    --------
    - count: array
        ann array of size (3 x nb proteins), index 0 corresponding to 
        the number of hydrophile, index 1 to the number of hydrophobic
        and index 2 to the number of gaps inside a hydrophobic cluster 
        at each columns of the MSA
    """
    count = np.zeros((2, size))
    nb_prot = len(clusters)
    for prot in clusters:
        for cl in clusters[prot]:
            start = cl[0]
            for i in range(len(cl[2])):
                #count[cl[2][i]] += 1
                if cl[2][i] == 0:
                    count[0, start+i] += 1
                elif cl[2][i] == 1:
                    count[1, start+i] += 1
                elif cl[2][i] == -1:
                    continue
                    #count[2, start+i] += 1
                else:
                    raise ValueError("The impossible happened, where does this value come from: {}! (protein {})".format(cl[2][i], prot))
    count /= nb_prot
    return count

def filter_clusters(cluster_ali, conservation, cutoff):
    """ filter columns based on conservation, mark them as 
    True if kept else False    
    """
    filtered_clusters = dict()
    for prot in cluster_ali:
        clusters = list()
        for cl in cluster_ali[prot]:
            new_cl = list()
            start = cl[0]
            for i in range(len(cl[2])):
                if np.any(conservation[:,start+i] > cutoff):
                    new_cl.append(True)
                else:
                    new_cl.append(False)
            clusters.append(new_cl)
        filtered_clusters[prot] = clusters
    return filtered_clusters

def select_clusters(cluster_positions, conservation, cutoff):
    proteins = list(cluster_positions.keys())
    # temporary matrix, TODO replace cluster_positions directly by matrix
    mat_positions = np.zeros((len(proteins), len(cluster_positions[proteins[0]])))
    for i, prot in enumerate(proteins):
        mat_positions[i] = cluster_positions[prot]

    seeds = list()
    seed = None
    for col in range(mat_positions.shape[1]):
        if mat_positions[:, col].sum() > 0 and np.any(conservation[:, col] > cutoff):
            # cluster found:and conserved
            # is it a H, P or G ?
            idx = np.argmax(conservation[:, col])
            if idx == 2:
                idx = "*"
            idx = "{}".format(idx)
            # inside a seed
            #else:
            #    idx = "*"
            if seed is not None:
                seed.append(idx)
            elif seed is None and idx == "1":
                seed = [idx]
                start = col    
        else:
            if seed is not None:
                seeds.append((start, start+len(seed), seed))
            seed = None
        #print(col, "".join(seed) if seed is not None else None, conservation[:, col])
    if seed is not None:
        seeds.append((start, start+len(seed), seed))
    return seeds

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="inputfasta", 
                        help="input alignment")
    parser.add_argument("-c", action="store", dest="conservation",
                        help="conservation cutoff", type=float)
    parser.add_argument("-o", action="store", dest="outputfile", 
                        help="output file")
    params = parser.parse_args()
    return params

def main():
    params = get_cmd()
    
    # read input sequences
    dfasta = read_multifasta(params.inputfasta)
    check_msa(dfasta)

    queries = [q for q in dfasta.keys()]
    sequences = [str(dfasta[q].seq) for q in queries]
    seq_no_gap = [seq.replace("-", "").upper() for seq in sequences]

    # HCA
    hca = HCA(seq=seq_no_gap, querynames=queries)
    clusters = hca.get_clusters()

    # map clusters to sequences with gap
    clusters_idx, clusters_ali = map_clusters_ali(dfasta, clusters)
    
    # filter columns by conservation
    conservation = cluster_conservation(clusters_ali, len(sequences[0]))
    seeds = select_clusters(clusters_idx, conservation, params.conservation)
    
    with open(params.outputfile, "w") as outf:        
        for start, stop, seed in seeds:
            outf.write("{}\t{}\t{}\n".format(start+1, stop, "".join(seed)))
    
    sys.exit(0)

if __name__ == "__main__":
    main()
