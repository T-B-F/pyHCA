![Python 3.5](https://img.shields.io/badge/Python-3.5-blue.svg)

pyHCA and hcatk (HCA toolkit) are a python library and executable for the Hydrophobic Cluster Analysis of protein sequences.
pyHCA implements various class and function for handling sequences and analyses them
hcatk provides off the hands set of program to perform hydrophobic cluster analyses.

Requires
========

- Biopython >= 1.65
- ete3
- python3   >= 3.5

Installation
============

A quick install can be perform using:

    pip3 install .
 
However, the ete3 can be difficult to install as some features requires PyQt4 and sip.
Please refer to the official ete3 insallation guidelines [http://etetoolkit.org/download/] (http://etetoolkit.org/download/) for any support.
On Mac OS X, you will also need to install XQuartz to use ete3, please refer to [XQuartz documentation] (http://www.xquartz.org/).

We recommend you to work on a conda virtual environment to properly build the non Python extention of the ete3 package and afterward install pyHCA in this new environment.

Example
*******

download and install conda from miniconda website using the correct installer (64 bits / 32 bits, MacOSC / Linux):

on MacOSX

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh

on Linux (64 bits installer)

    https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

create a virtual environment and switch to the environment

    conda create -n test_pyHCA python=3.5 pip
    source activate test_pyHCA
    
and install pyqt4 before running pyHCA installer

    conda install -n test_pyHCA pyqt=4.11.4
    cd <path to pyHCA directory>
    pip install .
    


Usage
=====

segment
-------

Use the composition in hydrophobic cluster of a sequence to detect domains.

    $ hcatk segment -h

    usage: hcatk segment [-h] -i INPUTF -o OUTPUTF [-v]
                          [-m {cluster,domain}]
                          [-t {aminoacid,nucleotide}]

Arguments:
**********

    -h, --help            show help message and exit
    
required arguments:

    -i INPUTF             an amino-acid sequence files in fasta format
    -o OUTPUTF            the output file with annotation
    
optional arguments:

    -v                    keep temporary results
    -m {cluster,domain}   method to use, cluster: will report *the hydrophobic
                          clusters found in the sequence, domain: will delineate
                          domains based on the hydrophobic cluster profile of
                          the sequence
    -t {aminoacid,nucleotide}
                          the type of the biological sequences passed in the
                          input file

Example:
********

    $ hcatk segment -i data/orc1.fasta -o data/orc1.hca -m domain


Format:
*******

The output is formated in a fasta like style, with an header storing the protein name and size:

    >protein_name protein_size
    
followed by four columns:

    domain  524     527     nan
    domain  552     923     0.0032921164246364487
    cluster 1       2       11
    cluster 10      17      10001011
    cluster 23      23      1
    
The first column correspond to the hca element identified, either a domain or a cluster.
The second and third columns correspond to the start and stop (indexed from 1 to the sequence length) of the hca element.
The fourth column either corresponds to a p-value if the element is a hca domain or to the hydrophobic cluster in binary mode if the element is a cluster.
The p-value of the domain is computed against a reference distribution made of intrinsically disordered proteins and describe the "degree of foldability associated to the hca domain element.
A nan value is returned if the domain has less than 30 residues

    
draw
----

Draw the HCA diagram of a sequence.
Optionnaly, can display the domain annotation of a sequence if provied.

    $ hcatk draw -h

    
    usage: hcatk draw [-h] -i FASTAFILE [-w WINDOW] [-d DOMAIN] [-f {pfam,seghca}]
                  [--color-msa {rainbow,identity}] -o OUTPUTFILE

Arguments:
**********

    -h, --help            show help message and exit
    
required arguments:
    
    -i FASTAFILE          the fasta file
    -o OUTPUTFILE         output file in svg format

optional arguments:

    -w WINDOW             sequence len before breaking the sequence to the next
                          plot (-1 the whole sequence are used, minimum size is
                          80)
    -d DOMAIN             [optionnal] provide domain annoation
    -f {pfam,seghca}      the domain file format
    --cons-msa {aa,hca}   method to use to compare sequences (aa, hca)

Example:
********
    
    $ hcatk draw -i data/PF00533_sub.txt -o data/PF00533_sub.svg --cons-msa aa
    $ inkscape data/PF00533_sub.svg # external svg viewer


tremolo
-------

Use Tremolo-HCA to find remote homologous proteins with domain context.


    $ hcatk tremolo -h

    usage: hcatk [-h] -f INPUTFASTA [-d DOMAINS [DOMAINS ...]] 
                 [-w WORKDIR] [-E EVALUE] [-o OUTPUT] 
                 [--hhblits-params HHBLITSPARAMS] [--hhblits-db HHBLITSDB]


Arguments:
**********

    -h, --help            show help message and exit

required arguments:

    -i INPUTFASTA         input fasta file
    -o OUTPUT             output file
    -w WORKDIR            working directory

optional arguments:

    -d DOMAINS [DOMAINS ...] list of domain positions (start and stop inclusive
                             and separated by comma : -d 1,10 20,30 60,100. 
                             If not provided the search will be performed on 
                             each domain found after segmentation of the input 
                             sequence. To use the whole protein use -d all.
    --p2ipr P2IPR            path to the Interpro annotation of UniproKBt proteins,
                             gzip format supported.
    -E EVALUE                filter hhblits results by evalue
    --hhblits-params HHBLITSPARAMS 
                            parameters to pass to hhblits, between quotes
    --hhblits-db HHBLITSDB  path to the database to use with hhblits

Example:
********

    $ hcatk tremolo -i data/orc1.fasta --p2ipr data/protein2ipr.dat.gz -E 0.001 --hhblits-db hhsuite/uniprot20_2016_02/uniprot20_2016_02 -o data/orc1_tremolo.txt -w tremolo_tmp


domOnTree
---------

Vizualise Tremolo-HCA results on a taxonomic tree with protein domain arrangement information


    $ hcatk domOnTree -h

    usage: hcatk [-h] [-i TREMOLORES] [-t TREEFILE] [-s PROT2SPECIES]
                 [-n NCBITAXID [NCBITAXID ...]] [-o OUTPUT]


Arguments:
**********

    -h, --help            show help message and exit
    
required arguments:

    -i TREMOLORES         tremolo results with domain matchs
    -o OUTPUT             phylogenetic tree with tremolo hits
    
optional arguments:


    -t TREEFILE           phylogenetic tree with node as ncbi taxonomic ids
    -s PROT2SPECIES       file with prot to species informations
    -n NCBITAXID [NCBITAXID ...]
                          list of node for which leaves will be merged (internal
                          node need to be in tree)


Example:
********

    $ hcatk domOnTree -i data/orc1_tremolo.txt -o data/orc1_tremolo.pdf
                          

Additional ressources
---------------------

The interpo domain annoation can be downloaded at:
wget ftp://ftp.ebi.ac.uk/pub/databases/interpro/protein2ipr.dat.gz

HHblits of the HH-suite package can be downlad at (v3 or higher):
git clone git@github.com:soedinglab/hh-suite.git

And the uniprot hhblits compatible database at:
http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
