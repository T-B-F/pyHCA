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
from Bio.Alphabet import IUPAC
import numpy as np
from .annotateHCA import _annotation_aminoacids
from .drawHCA import make_svg, getSVGheader, colorize_positions
from .seq_util import compute_conserved_positions

class HCA(object):
    """ HCA class provides an API interface to all the standalone programs
    """
    def __init__(self, querynames=None, seq=None, seqfile=None, file_format="fasta"):
        """ create an HCA instance, accept only one sequence at a time
        
        Parameters:
        -----------
        seq: list of instances or instance
            instance can either be a string, a Bio.Seq object or a Bio.SeqRecord object
        querynames: list or string
            if seq is instance of string or Bio.Seq, provide a name to the protein sequences, if not used "query_<num_idx>" is used
        seqfile: string
            path to the file containing protein sequences
        file_format: string
            BioPython supported file format for seqfile
            
            
        Usage:
        ------
        >>> # instanciation of a single sequence
        >>> seq_str = "ATGYHVVLIVQEAGFHILLV"
        >>> hca = HCA(seq_str, querynames="my_query")
        >>>        
        >>> from Bio import Seq
        >>> seq_bio = Seq.Seq("ATGYHVVLIVQEAGFHILLV")
        >>> hca = HCA(seq_bio) # without specifying querynames, automatically set it up to query_1
        >>>
        >>> from Bio import SeqRecord
        >>> seq_rec = SeqRecord.SeqRecord(id="protein_1", seq="ATGYHVVLIVQEAGFHILLV")
        >>> hca = HCA(seq_rec) # SeqRecord's id attribute is used 
        >>>

        >>> # instanciation of a list of sequences
        >>> seq_str_list = ["ATGYHVVLIVQEAGFHILLV", "AGVVLATGYHHILLVFHILLV"]
        >>> hca = HCA(seq_str_lst, querynames=["my_query_1", "my_query_2"])
        >>>
        """
        
        if is_seq_none(seq) and seqfile == None:
            raise ValueError("Error, you need to provide a sequence (or a list of sequences), or a valid path to a file to instantiate a HCA object")
        
        if not is_seq_none(seq) and seqfile != None:
            raise ValueError("Error, you need to provide either a sequence (or a list of sequences), or a valid path to a file to instantiate a HCA object")
        
        self.sequences = list()
        self.msa_seq = list()
        self.is_msa = False
        self.querynames = list()
        self.__domains = dict()
        self.__clusters = dict()
        self.__scores = dict()
        self._number_of_sequences = 0
        
        # attributes to keep in memory if computation were previously done
        self._segments_done = False
        self._segments_done_with_t = -1
        self._tremolo_done = False
        
        if not is_seq_none(seq):
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
                    
        # check if sequence is a MSA sequence
        is_msa, msa_seq, sequences = check_if_msa(self.querynames, self.sequences)
        self.msa_seq = msa_seq[:]
        self.is_msa = is_msa
        if is_msa:
            self.sequences = sequences[:]
        
    ### SEG-HCA
    @property
    def domains(self):
        return self.__domains
    
    @property
    def clusters(self):
        return self.__clusters
    
    @property
    def scores(self):
        return self.__scores
    
    @domains.setter
    def domains(self, domains):
        self.__domains = domains
    
    @clusters.setter
    def clusters(self, clusters):
        self.__clusters = clusters
    
    @scores.setter
    def scores(self, scores):
        self.__scores = scores
        
    def segments(self, t=0.1, verbose=False):
        """ run the segmentation of protein sequences into HCA domains, store domain and cluster positions
        """
        t1 = time.time()
        if not self._segments_done or t != self._segments_done_with_t:
            if t != self._segments_done_with_t and verbose:
                print("Running segmentation with a different t value ({} -> {})".format(self._segments_done_with_t, t))
            for i in range(len(self.sequences)):
                seq = self.sequences[i]
                prot = self.querynames[i]
                annotations = _annotation_aminoacids(seq, t=t, method="domain", verbose=verbose)
                self.domains[prot] = annotations["domain"]
                self.clusters[prot] = annotations["cluster"]
                self.scores[prot] = annotations["scores"]
                self._segments_done = True
                self._segments_done_with_t = t
        elif verbose:
            print("Warning, segmentations already performed return precomputed results, "
                    "used a different t value ({}) to get different results".format(self._segments_done_with_t), file=sys.stderr)
        t2 = time.time()
        if verbose:
            print("Segmentation done in {}".format(t2-t1))
    
    def get_domains(self, prot=None):
        """ function wrapper to return HCA domain annotation.
        If only one sequence was provided return a list of domains.
        If multiple sequences were provided return a dictionary with
        protein quernames as keys and the list of domains as values.
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.domains:
                raise KeyError("Error, unable to find proteins '{}' in domain results".format(prot))
            return self.domains[prot]
        elif self._number_of_sequences == 1:
            return self.domains[self.querynames[0]]
        else:
            return [self.domains[prot] for prot in self.querynames]
    
    def get_clusters(self, prot=None):
        """ function wrapper to return HCA cluster positions. 
        If only one sequence was provided return a list of clusters.
        If multiple sequences were provided return a dictionary with
        as protein querynames as keys and the list of clusters as values.
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.clusters:
                raise KeyError("Error, unable to find proteins '{}' in cluster results".format(prot))
            return self.clusters[prot]
        elif self._number_of_sequences == 1:
            return self.clusters[self.querynames[0]]
        else:
            return [self.clusters[prot] for prot in self.querynames]
        
    
    def get_scores(self, prot=None):
        """ function wrapper to return HCA scores of each domains. 
        If only one sequence was provided return a list of scores.
        If multiple sequences were provided return a dictionary with
        as protein querynames as keys and the list of scores as values.
        """
        if not self._segments_done:
            self.segments()
        if prot != None:
            if prot not in self.scores:
                raise KeyError("Error, unable to find proteins '{}' in scores results".format(prot))
            return self.scores[prot]
        elif self._number_of_sequences == 1:
            return self.scores[self.querynames[0]]
        else:
            return [self.scores[prot] for prot in self.querynames]
    
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
        if self.is_msa:
            msa_conserved = compute_conserved_positions(dict(zip(self.querynames, self.sequences)), dict(zip(self.querynames, self.msa_seq)))
        
        self.all_svg = dict()
        max_aa = 0
        cnt = 0
        # create hca plot for each sequence
        for i in range(len(self.querynames)):
            prot = self.querynames[i]
            prot_seq = self.sequences[i]
            # read domain annotation if provided
            annotation = {"domains": list()}
            if prot in external_annotation and "domains" in external_annotation[prot]:
                for start, stop, dom in external_annotation[prot]:
                    annotation["domains"].append((start, stop, dom, "!", None))
            if show_hca_dom:
                for dom in self.domains[prot]:
                    start = dom.start
                    stop = dom.stop
                    annotation["domains"].append((start, stop, "domain", "!", None))
            # make svg
            
            if self.is_msa:
                annotation["positions"] = colorize_positions(self.msa_seq[i], prot_seq, msa_conserved[prot], method="rainbow")
            cur_svg, nbaa = make_svg(prot, prot_seq, annotation, cnt)
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
    
def is_seq_none(seq):
    """ check if seq argument is None
    """
    # Bio.SeqRecord.SeqRecord does not support direct comparison
    if isinstance(seq, Bio.SeqRecord.SeqRecord):
        return False
    return seq == None
    
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

def check_if_msa(querynames, sequences):
    """ check if provided sequences are from an MSA
    if yes, transform msa sequences to ungapped sequences for HCA analysis
    """
    msa_length, msa_seq, new_sequences = list(), list(), list()
    is_msa = False
    prot_alphabet = set(IUPAC.protein.letters)
    gaps = set(["-", "."])
    for i, seq in enumerate(sequences):
        new_seq = ""
        for j, c in enumerate(seq):
            if c in prot_alphabet:
                new_seq += c
            else:
                if c in gaps:
                    is_msa = True
                else:
                    print("Invalid amino acids ({}, {}) in protein {}, replaced by X".format(j, c, querynames[i]))
        new_sequences.append(new_seq)
        msa_seq.append(seq)
        msa_length.append(len(seq))
        
    if is_msa:
        # chec identical sequence lengths
        if len(set(msa_length)) != 1:
            raise ValueError("Error, MSA characters detected but sequences have different lengths")
    
    return is_msa, msa_seq, new_sequences
