#!/usr/bin/env python
""" The classHCA is an ensemble of class to describe hydrophobic clusters and
domain created through hydrophobic clusters.
"""

__author__ = "Tristan Bitard-Feildel"
__licence__= "MIT"
__version__ = 0.1
__email__ = "tristan.bitard-feildel [you know what] impmc.upmc.fr"
__institute__ = "IMPMC, UPMC"

from Bio import Seq 

class HCA(object):
    
    def __init__(seq=None, seqtype=None, seqfile=None, file_format="fasta"):
        if seq == None and seqfile == None:
            raise ValueError("Error, you need to provide a sequence or a valid path to a file to instantiate a HCA object")
        
        if seq != None and seqfile != None:
            raise ValueError("Error, you need to provide either a sequence or a valid path to a file to instantiate a HCA object")
        
        if seq != None:
            # it's a sequence
            if isinstance(seq, str):
                self.seq = seq
            if isinstance(seq, Seq.Seq):
                self.seq = str(seq)
            elif isinstance(seq, Bio.SeqRecord.SeqRecord):
                self.seq = str(seq.seq)
            else:
                raise TypeError("Error, the 'seq' argument must be either a string a Bio.Seq instance or a Bio.SeqRecord instance")
        
        else:
            # it's a file
            if not os.path.isfile(seqfile):
                raise ValueErro("Error unable to find file {}".format(seqfile))
            with open(seqfile) as inf:
                record = SeqIO.read(inf, file_format)
            self.seq = str(record.seq)
            
            self._segments_done = False
            self._tremolo_done = False
            
        def segments(self):
            self._segments_done = True
            pass
        
        def get_domains(self):
            if not self._segments_done:
                self.segments()
            return self.domains
        
        def get_clusters(self):
            if not self._segments_done:
                self.segments()
            return self.clusters
        
        def tremolo(self):
            pass
        
        def draw(self):
            pass
        
        
        