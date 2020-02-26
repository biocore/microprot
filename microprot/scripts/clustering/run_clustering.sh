#!/bin/bash

help=" Script clusters faa file using MMseqs in three steps:
 - creates database
 - uses linclust for initial clustering
 - uses clust for final clustering

 Output:
 - fasta file from linclust
 - fasta file from clust
 - fasta file from clust with cluster sizes added to headers
 - tsv file for statistical analysis with the structure:
   sequence_length\tcluster_size

 Requirements:
 - Python 3.x
 - MMseqs2 added to PATH

 Usage:
   [-h] show help
   [-D] create DB
   [-L] run linclust
   [-C] run clust
   [-i] full path to the input file
   [-m] --min-seq-id    for linclust
   [-c] --cov-mode      for linclust
   [-k] --kmer-per-seq  for linclust
   [-m] --min-seq-id    for clust

 Example:
   ./run_clustering.sh -DLC -i /full/path/to/database/db.faa -m 0.7 -c 1 -k 80 -M 0.3
"

timestamp() {
    date --rfc-3339=seconds
}

export -f timestamp

OPTIND=1
while getopts "h?DLCi:m:c:k:M:" opt; do
    case "$opt" in
    h|\?)
        echo "${help}"
        exit 0
        ;;
    D)  CREATE_DB=1
        ;;
    L)  RUN_LINCLUST=1
        ;;
    C)  RUN_CLUST=1
        ;;
    i)  export INPUT_PATH=$OPTARG
        ;;
    m)  export MINSEQID_LINCLUST=$OPTARG
        ;;
    c)  export COVMODE_LINCLUST=$OPTARG
        ;;
    k)  export KMERS_LINCLUST=$OPTARG
        ;;
    M)  export MINSEQID_CLUST=$OPTARG
        ;;
    esac
done

export DB_INPUT_FILE=$(basename "${INPUT_PATH}")
export MAINDIR=$(dirname "${INPUT_PATH}")

# File identifiers
export IDENTIFIER_LINCLUST="${MINSEQID_LINCLUST}_\
${COVMODE_LINCLUST}_${KMERS_LINCLUST}"
export IDENTIFIER="${IDENTIFIER_LINCLUST}__${MINSEQID_CLUST}"

# Create DB for the input file
if [[ "${CREATE_DB}" = 1 ]]; then
    ./create_db.sh
fi

# Run linclust
if [[ "${RUN_LINCLUST}" = 1 ]]; then
    ./linclust.sh
fi

# Run clust for linclust results from above
if [[ "${RUN_CLUST}" = 1 ]]; then
    ./clust.sh
fi


