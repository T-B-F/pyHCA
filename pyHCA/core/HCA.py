#!/usr/bin/env python
""" The classHCA is an ensemble of class to describe hydrophobic clusters and
domain created through hydrophobic clusters.
"""

__author__ = "Tristan Bitard-Feildel"
__licence__= "MIT"
__version__ = 0.1
__email__ = "tristan.bitard-feildel [you know what] impmc.upmc.fr"
__institute__ = "IMPMC, UPMC"

import time, os
import Bio
from Bio import Seq 
from Bio import SeqIO
from .annotateHCA import _annotation_aminoacids
from .drawHCA import make_svg, getSVGheader

class HCA(object):
    """ HCA class provides an API interface to all the standalone programs
    """
    def __init__(self, querynames=None, seq=None, seqtype=None, seqfile=None, file_format="fasta"):
        """ create an HCA instance, accept only one sequence at a time
        """
        
        if seq == None and seqfile == None:
            raise ValueError("Error, you need to provide a sequence (or a list of sequences), or a valid path to a file to instantiate a HCA object")
        
        if seq != None and seqfile != None:
            raise ValueError("Error, you need to provide either a sequence (or a list of sequences), or a valid path to a file to instantiate a HCA object")
        
        self.sequences = list()
        self.querynames = list()
        self.__domains = dict()
        self.__clusters = dict()
        
        self._number_of_sequences = 0
        self._segments_done = False
        self._tremolo_done = False
        self._segments_done_with_t = -1
        
        if seq != None:
            # list of sequences or single sequence
            qnames = list()
            if isinstance(seq, list):
                for s in seq:
                    new_seq, qname = check_seq_type(s)
                    self.sequences.append(new_seq)
                    qnames.append(qname)
                    self._number_of_sequences += 1
            else:
                seq, qname = check_seq_type(seq)
                self.sequences.append(seq)
                qnames = [qname]
                self._number_of_sequences += 1
            
            # adapt queryname to list type
            if querynames == None:
                all_qnames = list(set(qnames))
                if len(all_qnames) == 1 and all_qnames[0] == "query":
                    for idx in range(len(self.sequences)):
                        name = "query_{}".format(idx+1)
                        self.querynames.append(name)
                        # initialize containers for domains and clusters
                        if name not in self.__domains:
                            self.__domains[name] = list()
                            self.__clusters[name] = list()
                        else:
                            raise RuntimeError("Error, multiple proteins found with the same name {}".format(name))
                            
                else:
                    for name in qnames:
                        self.querynames.append(name)
                        # initialize containers for domains and clusters
                        if name not in self.__domains:
                            self.__domains[name] = list()
                            self.__clusters[name] = list()
                        else:
                            raise RuntimeError("Error, multiple proteins found with the same name {}".format(name))
            elif isinstance(querynames, list):
                if len(querynames) != len(self.sequences):
                    raise ValueError("Error, queryname argument must be a list of same size as seq argument (one name per seq")
                else:
                    for name in querynames:
                        self.querynames.append(name)
                        # initialize containers for domains and clusters
                        if name not in self.__domains:
                            self.__domains[name] = list()
                            self.__clusters[name] = list()
                        else:
                            raise RuntimeError("Error, multiple proteins found with the same name {}".format(name))
            elif isinstance(querynames, str):
                self.querynames.append(querynames)
            else:
                raise ValueError("Error, querynames should be a list of names or a string")
        else:
            # it's a file
            if not os.path.isfile(seqfile):
                raise ValueErro("Error unable to find file {}".format(seqfile))
            with open(seqfile) as inf:
                for record in SeqIO.parse(inf, file_format):
                    self.sequences.append(str(record.seq))
                    self.querynames.append(record.id)
                    self._number_of_sequences += 1
                    # initialize containers for domains and clusters
                    if record.id not in self.__domains:
                        self.__domains[record.id ] = list()
                        self.__clusters[record.id ] = list()
                    else:
                        raise RuntimeError("Error, multiple proteins found with the same name {}".format(record.id ))
        
    ### SEG-HCA
    @property
    def domains(self):
        return self.__domains
    
    @property
    def clusters(self):
        return self.__clusters
        
    @domains.setter
    def domains(self, domains):
        self.__domains = domains
    
    @clusters.setter
    def clusters(self, clusters):
        self.__clusters = clusters
        
    def segments(self, t=0.1, verbose=False):
        """ run the segmentation in domains, store domain and cluster positions
        """
        t1 = time.time()
        if not self._segments_done or t != self._segments_done_with_t:
            for i in range(len(self.sequences)):
                seq = self.sequences[i]
                prot = self.querynames[i]
                annotations = _annotation_aminoacids(seq, t=t, method="domain", verbose=verbose)
                self.domains[prot] = annotations["domain"]
                self.clusters[prot] = annotations["cluster"]
                self._segments_done = True
                self._segments_done_with_t = t
        elif verbose:
            print("Warning, segmentations already performed return precomputed results, "
                    "used a different t value ({}) to get different results".format(self._segments_done_with_t), file=sys.stderr)
        t2 = time.time()
        if verbose:
            print("Segmentation done in {}".format(t2-t1))
    
    def get_domains(self, prot=None):
        """ get domain annotations
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.domains:
                raise KeyError("Error, unable to find proteins '{}' in domain results".format(prot))
            return self.domains[prot]
        return [self.domains[prot] for prot in self.querynames]
    
    def get_clusters(self, prot=None):
        """ get cluster positions
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.clusters:
                raise KeyError("Error, unable to find proteins '{}' in cluster results".format(prot))
            return self.clusters[prot]
        return [self.clusters[prot] for prot in self.querynames]
    
    def save_annotation(self, output):
        """ save the seg-HCA annotation to a file
        """
        if not self._segments_done:
            raise RuntimeError("Error, no annotation to save, you must perform an HCA segmentation with the segments() method first")
        with open(output, "w") as outf:
            for i, prot in enumerate(self.querynames):
                outf.write(">{} {}\n".format(prot, len(self.sequences[i])))
                for domannot in self.domains[prot]:
                    outf.write("{}\n".format(str(domannot)))
                for clustannot in self.clusters[prot]:
                    outf.write("{}\n".format(str(clustannot)))
            
    ### DrAW-HCA
    def draw(self, external_annotation=dict(), show_hca_dom=False, outputfile=None):
        """ draw a HCA plot in svg of each sequence
        """
        self.all_svg = dict()
        max_aa = 0
        cnt = 0
        # create hca plot for each sequence
        for i in range(len(self.querynames)):
            prot = self.querynames[i]
            prev_seq = self.sequences[i]
            # read domain annotation if provided
            annotation = list()
            if prot in external_annotation:
                for start, stop, dom in external_annotation[prot]:
                    annotation.append((start, stop, dom, "!", None))
            if show_hca_dom:
                for dom in self.domains[prot]:
                    start = dom.start
                    stop = dom.stop
                    annotation.append((start, stop, "domain", "!", None))
            # make svg
            cur_svg, nbaa = make_svg(prot, prev_seq, annotation, cnt)
            self.all_svg[prot] = cur_svg
            if nbaa > max_aa:
                max_aa = nbaa
            cnt += 1
        # write in outputfile if provided
        if outputfile != None:
            svgheader = getSVGheader(max_aa, (cnt+1)*230)
            with open(outputfile, "w") as fout:
                fout.write(svgheader)
                for i in range(len(self.querynames)):
                    prot = self.querynames[i]
                    fout.write(self.all_svg[prot])
                fout.write("</svg>")
        # return the svg dictionary
        return self.all_svg
    
def check_seq_type(seq):
    """ check that sequence is of correct type
    """
    queryname = "query"
    if isinstance(seq, str):
        seq = seq
    elif isinstance(seq, Seq.Seq):
        seq = str(seq)
    elif isinstance(seq, Bio.SeqRecord.SeqRecord):
        seq = str(seq.seq)
        queryname = seq.id
    else:
        raise TypeError("Error, the 'seq' argument must be either a string a Bio.Seq instance or a Bio.SeqRecord instance")
    return seq, queryname
