import click


@click.command()
@click.option('--occ_file', '-c', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Path to linclust occ tsv file.')
@click.option('--clust_file', '-f', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Path to clust tsv file.')
@click.option('--out_file', '-o', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Path to file where output will be saved.')
def add_occ(occ_file, clust_file, out_file):
    """
    Add third column with cluster sizes (here: occurrences)
    to the clust output tsv file based on the linclust occ
    tsv file and save result to clust app tsv file.

    Exemplary input:

        linclust occ tsv file (occ_file):

        GUT_GENOME152271_02644 2
        GUT_GENOME250812_01207 1
        (...)

        clust output tsv file (clust_file):

        GUT_GENOME152271_03641 GUT_GENOME250812_01207
        GUT_GENOME224440_07043 GUT_GENOME152271_02644
        (...)

    Exemplary output (out_file):

        GUT_GENOME152271_03641 GUT_GENOME250812_01207 1
        GUT_GENOME224440_07043 GUT_GENOME152271_02644 2
        (...)
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
    add_occ()
