#!/bin/bash
#PBS -N microprot_test
#PBS -m a
#PBS -t 1-10%10
#PBS -l nodes=1:ppn=4
#PBS -l mem=50gb
#PBS -l walltime=100:00:00
#PBS -o /home/tkosciolek/cluster/microprot_logs/${PBS_JOBNAME}_${PBS_JOBID}.o
#PBS -e /home/tkosciolek/cluster/microprot_logs/${PBS_JOBNAME}_${PBS_JOBID}.e

set -e
uname -a

source activate microprot

# do I want to copy source FASTA file too?

workdir="/projects/microprot/results/"


JOBNUM=`echo ${PBS_JOBID} | cut -d '[' -f1`

tmpdir=/localscratch/microprot/${PBS_JOBNAME}_${JOBNUM}
rm -rf $tmpdir
mkdir -p ${tmpdir}
cd ${tmpdir}

ln -s ${workdir}/config.yml ${tmpdir}
date --rfc-3339=seconds

snakemake -s /projects/microprot/microprot/snakemake/Snakefile \
          --configfile ${tmpdir}/config.yml \
          --config seq_no=${PBS_ARRAYID} MICROPROT_TEMP=${tmpdir} \
	      --restart-times 1 --keep-going --latency-wait 5

date --rfc-3339=seconds
rm -r $tmpdir
