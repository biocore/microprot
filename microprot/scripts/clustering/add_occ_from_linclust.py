"""
Run calc_occ_linclust.py first!

Script adds third column with cluster sizes (here:
occurrences) to the clust output tsv file based on the
linclust occ tsv file and saves result to clust app
tsv file.

Exemplary input:

    linclust occ tsv file:

    GUT_GENOME152271_02644 2
    GUT_GENOME250812_01207 1
    (...)

    clust output tsv file:

    GUT_GENOME152271_03641 GUT_GENOME250812_01207
    GUT_GENOME224440_07043 GUT_GENOME152271_02644
    (...)

Exemplary output:

    GUT_GENOME152271_03641 GUT_GENOME250812_01207 1
    GUT_GENOME224440_07043 GUT_GENOME152271_02644 2
    (...)
"""

import sys


def add_occ(occ_file, clust_file, out_file):
    """
    Parameters
    ----------
    occ_file : path
        path to linclust occ tsv file
    clust_file : path
        path to clust tsv file
    out_file : path
        path to file where output will be saved
    """
    print("Reading lines...")
    with open(occ_file, "r") as datafile:
        data_occ = datafile.readlines()

    print("Creating dict...")
    dd = {}
    for line in data_occ:
        key, val = line.split()
        dd[key] = val

    print("Joining...")
    with open (clust_file, "r") as datafile:
        data_to_join = datafile.readlines()

    occ_sum = 0  # sum of occurrences
    res = []  # list of representative sequences
    for i, line in enumerate(data_to_join):
        seq_rep, seq_to_find = line.split()
        occ = dd[seq_to_find]
        occ_sum += int(occ)
        res.append(seq_rep + " " + seq_to_find + " " + occ + "\n")

    print("Writing results...")
    with open (out_file, "w") as resfile:
        resfile.write("".join([el for el in res]))

    print("Number of occurrences: {}\n".format(occ_sum))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Wrong number of arguments!")
    else:
        occ_file = sys.argv[1]
        clust_file = sys.argv[2]
        out_file = sys.argv[3]
        add_occ(occ_file, clust_file, out_file)
