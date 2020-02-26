#!/bin/bash

timestamp
START=$(date +%s.%N)

INPDIR="${MAINDIR}"
DBDIR="${MAINDIR}/output_${IDENTIFIER_LINCLUST}"
OUTDIR="${MAINDIR}/output_${IDENTIFIER}"
TMPDIR="${MAINDIR}/tmp_${IDENTIFIER}"
DB_NAME="${DB_INPUT_FILE%.*}_${MINSEQID_LINCLUST}"
SEQUENCE_DB="${DBDIR}/${DB_NAME}"

mkdir -p "${TMPDIR}"
mkdir -p "${OUTDIR}"

#============================
echo -e "\nPREFILTER...\n"
#============================
INPUT="${SEQUENCE_DB}"
TMPPATH="${TMPDIR}"
OUTPUT="${TMPPATH}/prefilter"
PARAMS="--max-seqs 300 -c 0.8 --comp-bias-corr 1 -s 6"
mmseqs prefilter "${INPUT}" "${INPUT}" "${OUTPUT}" ${PARAMS}
timestamp

#============================
echo -e "\nALIGN...\n"
#============================
PREFILTER_RESULT="${OUTPUT}"
OUTPUT="${TMPPATH}/align"
PARAMS="-c 0.8 --alignment-mode 3 --min-seq-id ${MINSEQID_CLUST} \
 --comp-bias-corr 1 -e 0.001 --max-seq-len 32768 --max-rejected 2147483647"
mmseqs align "${INPUT}" "${INPUT}" "${PREFILTER_RESULT}" "${OUTPUT}" ${PARAMS}
timestamp

#============================
echo -e "\nCLUST...\n"
#============================
ALIGNMENT_RESULT="${OUTPUT}"
OUTPUT="${OUTDIR}/${DB_NAME}_${MINSEQID_CLUST}_clu"
PARAMS="--cluster-mode 0"
mmseqs clust "${SEQUENCE_DB}" "${ALIGNMENT_RESULT}" "${OUTPUT}" ${PARAMS}
mmseqs createtsv "${SEQUENCE_DB}" "${SEQUENCE_DB}" \
"${OUTPUT}" "${OUTPUT}.tsv"
timpestamp

#===============================
echo -e "\nCREATE OCC FILE...\n"
#===============================
PATH_LIN="${DBDIR}/${DB_NAME}"
python add_occ_from_linclust.py -c "${PATH_LIN}_clu_occ.tsv" \
-f "${OUTPUT}.tsv" -o "${OUTPUT}_app.tsv"

# Calculate cluster sizes (here: occurrences) from the
# output app tsv file and saves results to occ tsv file

awk 'BEGIN { occ = 0 }
{   curr = $1
    if(prev != curr && NR > 1){
        print prev "\t" occ
        occ = $3
    }
    else
        occ += $3
    prev = $1
}
END { print curr "\t" occ }' "${OUTPUT}_app.tsv" > "${OUTPUT}_occ.tsv"
timestamp

#=================================================
echo -e "\nEXTRACT REPRESENTATIVE SEQUENCES...\n"
#=================================================
ORDER_FILE="${OUTPUT}_order"
awk '{ print $1 }' "${OUTPUT}.index" > "${ORDER_FILE}"
OUTPUT="${OUTDIR}/${DB_NAME}_${MINSEQID_CLUST}"
mmseqs createsubdb "${ORDER_FILE}" "${INPUT}" "${OUTPUT}"
OUTPUT_FASTA="${INPDIR}/${DB_INPUT_FILE%.*}_${IDENTIFIER}.fasta"
mmseqs result2flat "${SEQUENCE_DB}" "${SEQUENCE_DB}" \
"${OUTPUT}" "${OUTPUT_FASTA}" "--use-fasta-header"

echo -e "\nAdd headers with cluster sizes to fasta file\n"

python add_occ_to_fasta.py -c "${OUTPUT}_clu_occ.tsv" -f "${OUTPUT_FASTA}" \
-o "${OUTPUT_FASTA%.*}_clulen.fasta" -s "${OUTPUT_FASTA%.*}_stat.tsv"
timestamp

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "Total time: " $DIFF " s"
echo "done."
