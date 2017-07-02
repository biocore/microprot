#!/bin/bash
#PBS -N microprot_array_test
#PBS -m a
#PBS -t 1-2%1
#PBS -l nodes=1:ppn=8
#PBS -l mem=50gb
#PBS -l walltime=3:00:00
#PBS -o /home/tkosciolek/cluster/microprot_logs/${PBS_JOBNAME}_${PBS_ARRAYID}.o${PBS_JOBID}
#PBS -e /home/tkosciolek/cluster/microprot_logs/${PBS_JOBNAME}_${PBS_ARRAYID}.e${PBS_JOBID}

set -e

source activate microprot
cd /projects/microprot/benchmarking/snakemake_minimal_test

# tmpdir=/localscratch/joseout/JOBNAME_${PBS_ARRAYID}
# mkdir -p ${tmpdir}

snakemake -s /projects/microprot/microprot/snakemake/Snakefile \
          --config seq_no=${PBS_ARRAYID} --configfile=`pwd`/config.yml

# rm -r $tmpdir
