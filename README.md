[![Coverage Status](https://coveralls.io/repos/github/biocore/microprot/badge.svg?branch=master)](https://coveralls.io/github/biocore/microprot?branch=master)
[![Build Status](https://travis-ci.org/biocore/microprot.svg?branch=master)](https://travis-ci.org/biocore/microprot)

# microprot
microProt is coded in Python 3.x

## Introduction
microProt clusters and annotates microbial metagenome sequences for the ultimate goal of predicting the 3-dimensional structure and function of these proteins.

## Install

## Requirements
Some of the tools and databases we're using were developed externally and cannot be automatically installed. We ask you to download them on your own, install and update appropriate paths in `paths.yml`

### dbs
* [PDB100](http://dunbrack.fccc.edu/Guoli/culledpdb_hh/pdbaanr.gz) from [PISCES](http://dunbrack.fccc.edu/PISCES.php) PDB culling server (updated weekly)
* [UniRef90](ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz) from EBI (the database is 10GB+; updated monthly)
* [Uniclust30](http://wwwuser.gwdg.de/~compbiol/uniclust/current_release/uniclust30_2016_09_hhsuite.tar.gz) (14GB+; updated every 3 months)
* [PfamA](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_30.0.tgz) for HH-suite (updated approx. every 6 months)

### tools
Tools requiring manual installation are listed and linked below:
* [HH-suite 3.0](https://github.com/soedinglab/hh-suite)
* [metaPSICOV](http://bioinfadmin.cs.ucl.ac.uk/downloads/MetaPSICOV/)
