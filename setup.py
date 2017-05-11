#!/usr/bin/env python3

import os, sys
from setuptools import setup, find_packages

__version__ = "0.1"


if sys.version_info < (3,4):
    print("python version must be at least 3.4")
    sys.exit(1)


# check if PyQt4 is here
try:
    import PyQt4
except ImportError:
    print("You must first install PyQt4 for ete3 to drawing functionnalities")
    sys.exit(1)

def readme():
    with open("README.md") as inf:
        return inf.read()

script_dir = os.path.join("pyHCA", "bin")
all_scripts = list()
for binf in ["hcatk"]: # add script names here
    all_scripts.append(os.path.join(script_dir, binf))

setup(name='pyHCA',
    version=__version__,
    author='Tristan Bitard-Feildel, Guillem Faure',
    author_email='tristan.bitard-feildel@impmc.upmc.fr',
    url='http://www.bornberglab.org/',
    description='Python library for HCA analysis',
    long_description=readme(),

    #classifiers=[
        #'Development Status :: 3 - Alpha',
        #'License :: OSI Approved :: MIT License',
        #'Programming Language :: Python :: 3',
        #'Topic :: Text Processing :: Linguistic',
        #],
    
    keywords=('protein, domain, HCA, '),
    license='MIT',
    scripts=all_scripts,
    packages=find_packages(exclude=[script_dir,]),
    include_package_data=True,
    install_requires=['biopython>=1.68',
                      'requests',
                      'six',
                      '
                      'ete3',
                      #'scikit-learn>=0.18.1',
                      #'numpy>=1.11.2'n
                      ],
    provides=['pyHCA({0:s})'.format(__version__)],
    zip_safe=False,
)

