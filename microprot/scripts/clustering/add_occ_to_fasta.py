"""
Script:
- adds cluster sizes (here: occurrences) to
  representative sequence headers in .fasta file
- creates tsv file for statistical analysis with
  the structure: sequence_length\tcluster_size
"""

import sys
import re


def add_occ(occ_file, fasta_file, out_fasta_file, out_stat_file):
    """
    Parameters
    ----------
    occ_file : path
        path to clust occ tsv file
    fasta_file : path
        path to fasta file
    out_fasta_file : path
        path to output fasta file with extended headers
    out_stat_file : path
        path to output file for statistical analysis
    """

    print("Reading cluster sizes...")
    with open(occ_file, "r") as datafile:
        data_occ = datafile.readlines()

    print("Reading fasta file...")
    with open(fasta_file, "r") as datafile:
        data_fasta = datafile.readlines()

    print("Creating dict...")
    dd = {}
    for line in data_occ:
        key, val = line.split()
        dd[key] = val

    print("Updating headers...")
    occ_sum = 0  # sum of occurrences
    res = []  # list of representative sequences
    stat = []  # list for statistical data
    for i in range(0, len(data_fasta), 2):
        header = data_fasta[i]
        seq = data_fasta[i+1]
        identifier = header.split()[0]
        # remove leading '>'
        identifier = re.sub("^>", "", identifier)
        occ = dd[identifier]
        occ_sum += int(occ)
        header = re.sub("\n$", " # " + occ + "\n", header)
        res.append(header)
        res.append(seq)
        # len must be calculated without "\n"
        seq = re.sub("\n$", "", seq)
        stat.append(str(len(seq)) + "\t" + occ + "\n")

    print("\nWriting results...")
    with open (out_fasta_file, "w") as resfile:
        resfile.write("".join([el for el in res]))
    with open (out_stat_file, "w") as resfile:
        resfile.write("".join([el for el in stat]))

    print("Number of sequences: {}\n".format(str(occ_sum)))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Wrong number of arguments!")
    else:
        occ_file = sys.argv[1]
        fasta_file = sys.argv[2]
        out_fasta_file = sys.argv[3]
        out_stat_file = sys.argv[4]
        add_occ(occ_file, fasta_file, out_fasta_file, out_stat_file)
