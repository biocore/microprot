#!/bin/bash

timestamp
START=$(date +%s.%N)

INPDIR="${MAINDIR}"
DBDIR="${MAINDIR}/output_db/"
OUTDIR="${MAINDIR}/output_${IDENTIFIER_LINCLUST}/"
TMPDIR="${MAINDIR}/tmp_${IDENTIFIER_LINCLUST}/"
DB_NAME="${DB_INPUT_FILE%.*}" # remove extension
SEQUENCE_DB="${DBDIR}/${DB_NAME}"

mkdir -p "${TMPDIR}"
mkdir -p "${OUTDIR}"

#==========================
echo -e "\nLINCLUST...\n"
#==========================
INPUT="${SEQUENCE_DB}"
OUTPUT="${OUTDIR}/${DB_NAME}_${MINSEQID_LINCLUST}_clu"
TMPPATH="${TMPDIR}"
PARAMS="--min-seq-id ${MINSEQID_LINCLUST} -c 0.8 \
--cov-mode ${COVMODE_LINCLUST} --comp-bias-corr 0 \
--cluster-mode 2 --kmer-per-seq ${KMERS_LINCLUST}"
mmseqs linclust "${INPUT}" "${OUTPUT}" "${TMPPATH}" ${PARAMS}
mmseqs createtsv "${SEQUENCE_DB}" "${SEQUENCE_DB}" \
"${OUTPUT}" "${OUTPUT}.tsv"
timestamp

#===============================
echo -e "\nCREATE OCC FILE...\n"
#===============================

# Calculate cluster sizes (here: occurrences) from the
# output tsv file and saves results to occ tsv file

awk 'BEGIN { occ = 0 }
{   curr = $1
    if(prev != curr && NR > 1){
        print prev "\t" occ
        occ = 1
    }
    else
        occ++
    prev = $1
}
END { print curr "\t" occ }' "${OUTPUT}.tsv" > "${OUTPUT}_occ.tsv"
timestamp

#=================================================
echo -e "\nEXTRACT REPRESENTATIVE SEQUENCES...\n"
#=================================================
ORDER_FILE="${OUTPUT}_order"
awk '{ print $1 }' "${OUTPUT}.index" > "${ORDER_FILE}"
OUTPUT="${OUTDIR}/${DB_NAME}_${MINSEQID_LINCLUST}"
mmseqs createsubdb "${ORDER_FILE}" "${INPUT}" "${OUTPUT}"
OUTPUT_FASTA="${INPDIR}/${DB_NAME}_${IDENTIFIER_LINCLUST}.fasta"
mmseqs result2flat "${SEQUENCE_DB}" "${SEQUENCE_DB}" \
"${OUTPUT}" "${OUTPUT_FASTA}" "--use-fasta-header"
timestamp

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "Total time: " $DIFF " s"
echo "done."
