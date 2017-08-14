[![Coverage Status](https://coveralls.io/repos/github/biocore/microprot/badge.svg?branch=master)](https://coveralls.io/github/biocore/microprot?branch=master)
[![Build Status](https://travis-ci.org/biocore/microprot.svg?branch=master)](https://travis-ci.org/biocore/microprot)

# microprot
microprot is coded in Python 3.x

## Introduction
microprot clusters and annotates microbial metagenome sequences for the ultimate goal of predicting the 3-dimensional structure and function of these proteins.

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

## Naming conventions

### Filenames
All filenames are in the form: `GenomeID`\_`GeneID`\_`ResiduesFrom`-`ResiduesTo`.  
For example, `CP003179.1_3319` means gene `3319` from genome `CP003179.1` (Sulfobacillus acidophilus DSM 10332), or `CP003179.1_3319_1-60` means amino acids 1 to 60 from that gene.

### File extensions

* a3m  
    An alignment file produced by HH-suite programs. It's a format similar to FASTA, but in sequence rows it contains additional information useful for the construction of HMMs (represented by [a-z]). A detailed description can be found in [HH-suite user guide](https://github.com/soedinglab/hh-suite/blob/master/hhsuite-userguide.pdf) (section 6.1).

* out  
    HH-suite output files reporting a list of hits for an input sequence, along with Probability, P-value, E-value and other parameters (hit list); as well as a set of pair-wise sequence alignments. A detailed description can be found in [HH-suite user guide](https://github.com/soedinglab/hh-suite/blob/master/hhsuite-userguide.pdf) (section 5).

* match  
    Internal `microprot` files showing which sub-sequence of the input sequence matched defined by `config.yml` criteria for any of `E-value`, `P-value`, `Prob` or `minimum sequence length` in the `.out` file. Multiple hits are possible. The file is reported in a FASTA format.

* non_match  
    All sub-sequences longer than the `minimum sequence length` that do not meet the criteria for `.match`. Internal `microprot` file.

### Example
Gene `CP00000.0_1` (`CP00000.0_1.fasta`) with 100 residues is run against HHsearch and it returns 2 outputs: `CP00000.0_1.out` and `CP00000.0_1.a3m`. Sequence split parameters are:
```
min_prob: 90.0
min_fragment_length: 10
```
and the hit list portion of `CP00000.0_1.out` is:
```
No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
 1 1ABC_A Uncharacterized protein  91.5   0.001   0.001   24.3   0.0   20   10-30    211-231 (260)
 2 1BCD_A Uncharacterized protein  90.3   0.001   0.001   26.4   0.0   55   33-88    28-83  (149)
 3 1CDE_A Uncharacterized protein  85.3     0.2   0.001   26.4   0.0   55   43-98    28-83  (149)
```

According to our criteria, hits 1 and 2 are matches (probability >= 90.0 and fragment length (from `Query_HMM`) >= 10).  
So `CP00000.0_1.match` file will contain sequences:
```
>CP00000.0_1_10-30
---------EXAMPLEEXAMPLEEXAMPL-----------------------------------------
------------------------------
>CP00000.0_1_33-88
---------------------------------EXAMPLEEXAMPLEEXAMPLEEXAMPLEEXAMPLEEX
AMPLEEXAMPLEEXAMPL------------
```
and `CP00000.0_1.non_match` will contain sequence:
```
>CP00000.0_1_89-100
----------------------------------------------------------------------
------------------EXAMPLEEXAMP
```
Sub-sequences `CP00000.0_1_1-9` and `CP00000.0_1_31-33` will be dropped from subsequent analyses, as they did not match `minimum fragment length` criteria.
