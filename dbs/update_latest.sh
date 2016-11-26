#!/bin/bash

# symlink all of the newest databases to the LATEST folder
### CHANGE THIS ###
NEWEST="/projects/microprot/dbs/Apr_2016"
###################

LATEST="/projects/microprot/dbs/LATEST"

# PDB70
rm $LATEST/pdb70/pdb*
for i in $NEWEST/pdb70/pdb*
do
    ln -s $i $LATEST/pdb70/
done

# UNIPROT20
rm $LATEST/uniprot20/uniprot20*
for i in $NEWEST/uniprot20/uniprot20*
do
    ln -s $i $LATEST/uniprot20/
done

# NR
rm $LATEST/nr/nr*
for i in $NEWEST/nr/nr*
do
    ln -s $i $LATEST/nr/
done
