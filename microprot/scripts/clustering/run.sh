#!/bin/bash

#============================================================
# Script clusters faa file using MMseqs in three steps:
# - creates database
# - uses linclust for initial clustering
# - uses clust for final clustering

# Usage:
# - change input parameters
# - uncomment appropriate commands
# - run script:
#                 ./run.sh
#============================================================

timestamp() {
    date --rfc-3339=seconds
}

export -f timestamp

#=================
# input parameters
#=================

# input file
export DB_INPUT_FILE="uhgp-100.faa"
# main directory where input file is stored
# and all output/tmp files will be saved
export MAINDIR="/klaster/work/pszczerbiak/uhgp-100_2020-01-31"
# path to MMseqs
MMSEQSPATH="/home/pszczerbiak/mmseqs/bin"
export PATH="$MMSEQSPATH:$PATH"
# linclust parameters
export MINSEQID_LINCLUST=0.7 # --min-seq-id
export COVMODE_LINCLUST=1    # --cov-mode
export KMERS_LINCLUST=80     # --kmer-per-seq
# clust parametes
export MINSEQID_CLUST=0.3    # --min-seq-id
# file identifiers
export IDENTIFIER_LINCLUST="${MINSEQID_LINCLUST}_\
${COVMODE_LINCLUST}_${KMERS_LINCLUST}"
export IDENTIFIER="${IDENTIFIER_LINCLUST}__${MINSEQID_CLUST}"

#=============================
# Create DB for the input file
#=============================

#./create_db.sh > "create_db.log" 2>&1

#=============
# Run linclust
#=============

#./run_linclust.sh > "run_linclust.log" 2>&1

#==========================================
# Run clust for linclust results from above
#==========================================

#./run_clust.sh  > "run_clust.log" 2>&1
