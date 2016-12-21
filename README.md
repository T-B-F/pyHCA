pyHCA and hcatk (HCA toolkit) are a python library and executable for the Hydrophobic Cluster Analysis of protein sequences.
pyHCA implements various class and function for handling sequences and analyses them
hcatk provides off the hands set of program to perform hydrophobic cluster analyses.

Requires
========
- Biopython >= 1.65
- python3   >= 3.4

Installation
============
pip3 install .

Usage
=====

annotate
--------

Use the composition in hydrophobic cluster of a sequence to detect domains.

$ hcatk annotate -h

usage: hcatk annotate [-h] -i INPUTF -o OUTPUTF [-v]
                      [-m {cluster,domain}]
                      [-t {aminoacid,nucleotide}]

Arguments:
~~~~~~~~~

  -h, --help            show this help message and exit
  -i INPUTF             an amino-acid sequence files in fasta format
  -o OUTPUTF            the output file with annotation
  -v                    keep temporary results
  -m {cluster,domain}
                        method to use, cluster: will report *the hydrophobic
                        clusters found in the sequence, domain: will delineate
                        domains based on the hydrophobic cluster profile of
                        the sequence
  -t {aminoacid,nucleotide}
                        the type of the biological sequences passed in the
                        input file

draw
----

Draw the HCA diagram of a sequence.
Optionnaly, can display the domain annotation of a sequence if provied.

$ hcatk draw -h

usage: hcatk draw [-h] -i FASTAFILE [-d DOMAIN] [-f {pfam,seghca}] -o SVGFILE

Arguments:
~~~~~~~~~

  -h, --help        show this help message and exit

required arguments:
  -i FASTAFILE      the fasta file
  -o SVGFILE        the svg file

optional arguments:
  -d DOMAIN         the domain file
  -f {pfam,seghca}  the domain file format


tremolo
-------

Use TremoloHCA to find remote homologous proteins with domain context.


$ hcatk tremolo -h

usage: hcatk [-h] -f INPUTFASTA [-d DOMAINS [DOMAINS ...]] [-w WORKDIR]
             [-E EVALUE] [-o OUTPUT] [--hhblits-params HHBLITSPARAMS]
             [--hhblits-db HHBLITSDB]

Arguments:
~~~~~~~~~~

  -h, --help            show this help message and exit

required arguments:
  -f INPUTFASTA         input fasta file
  -o OUTPUT             output file
  -w WORKDIR            working directory

optional arguments:
  -d DOMAINS [DOMAINS ...]
                        list of domain positions (start and stop inclusive and
                        separated by comma : -d 1,10 20,30 60,100. If not
                        provided the search will be performed on each domain
                        found after segmentation of the input sequence. To use
                        the whole protein use -d all.
  -a {CDD,Interpro}     defined annotation method to use (default=Interpro)
  --p2ipr P2IPR         path to the Interpro annotation of UniproKBt proteins.
                        If the argument is not specified and '-a Interpro' is
                        set, the annotation will be retrieve using web queries
                        of Biomart service which will be slower.
  -E EVALUE             filter hhblits results by evalue
  --hhblits-params HHBLITSPARAMS
                        parameters to pass to hhblits, between quotes
  --hhblits-db HHBLITSDB
                        path to the database to use with hhblits

The interpo domain annoation can be downloaded at:
wget ftp://ftp.ebi.ac.uk/pub/databases/interpro/protein2ipr.dat.gz
