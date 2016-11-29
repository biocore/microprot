# microprot
## Introduction
microProt clusters and annotates microbial metagenome sequences for the ultimate goal of predicting the 3 dimensional structure and function of these proteins.

## Install

## Requirements
Because we're using tools and databases developed externally, we ask you to download them on your own, install and update appropriate paths in `paths.yml`  

### dbs
* [PDB sequences](ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz); [PDB100 clusters](ftp://resources.rcsb.org/sequence/clusters/)
* [nr](ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz) (the database is 20GB+; also make sure to use BLAST `formatdb` on it, as specified in PSIPRED README file)
* [uniprot20](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/) (make sure to download the latest version of UniProt20)

### tools
* [HH-suite 3.0](https://github.com/soedinglab/hh-suite)
* [PSIPRED](http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/)
* [DISOPRED3](http://bioinfadmin.cs.ucl.ac.uk/downloads/DISOPRED/)
* [HMMer](http://hmmer.org/download.html)
* [metaPSICOV](http://bioinfadmin.cs.ucl.ac.uk/downloads/MetaPSICOV/)
