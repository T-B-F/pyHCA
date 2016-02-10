#!/usr/bin/env python

import os,sys, shlex, subprocess
from os.path import join
import platform
import distutils
from distutils import sysconfig
from setuptools import find_packages
from distutils.core import setup


__version__ = ''
with open('lib/pyHCA/__init__.py') as inp:
  for line in inp:
      if line.startswith('__version__'):
          exec(line.strip())
          break
          
PACKAGES = ['pyHCA']

PACKAGE_DIR = {}
for pkg in PACKAGES:
    PACKAGE_DIR[pkg] = join('lib', *pkg.split('.'))

SCRIPTS = ['scripts/hca_annotation']

if (platform.system() == 'Windows' or 
    len(sys.argv) > 1 and sys.argv[1] not in ('build', 'install')):
    for script in list(SCRIPTS):
        SCRIPTS.append(script + '.bat')

setup(
    name='pyHCA',
    version=__version__,
    author='Tristan Bitard-Feildel, Guillem Faure',
    author_email='t.bitard.feildel@uni-muenster.de',
    url='http://www.bornberglab.org/',
    description='Python library for HCA analysis',
    long_description="",
    license='MIT',
    keywords=('protein, domain, HCA, '),
    packages=PACKAGES,
    package_dir=PACKAGE_DIR,
    scripts=SCRIPTS,
    provides=['pyHCA({0:s})'.format(__version__)],
)




