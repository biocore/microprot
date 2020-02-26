#!/bin/bash

timestamp
START=$(date +%s.%N)

INPDIR="${MAINDIR}"
OUTDIR="${MAINDIR}/output_db/"
TMPDIR="${MAINDIR}/tmp_db/"
DB_NAME=${DB_INPUT_FILE}

mkdir -p "${TMPDIR}"
mkdir -p "${OUTDIR}"

#==========================
echo -e "\nCREATE DB...\n"
#==========================
INPUT="${INPDIR}/${DB_NAME}"
OUTPUT="${OUTDIR}/${DB_NAME%.*}" # remove extension
mmseqs createdb "${INPUT}" "${OUTPUT}"
timestamp

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo "Total time: " $DIFF " s"
