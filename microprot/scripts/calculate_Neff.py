import re
import math
import click
from scipy.cluster.hierarchy import linkage, fcluster
from skbio import io, Protein, TabularMSA, DistanceMatrix
from skbio.sequence.distance import hamming


"""
Calculate effective family size of a multiple sequence alignment with a given
percent identity cutoff.
"""


def parse_msa_file(infile):
    """Read sequences from a multiple sequence alignment (MSA) file.

    Parameters
    ----------
    infile : str
        file path to input MSA file in A3M format (like FASTA format, but
        lowercase letters will be dropped)

    Returns
    -------
    skbio TabularMSA
    """
    seqs = []
    for seq in io.read(infile, format='fasta'):
        seqs.append(Protein(re.sub('[a-z]', '', str(seq)),
                            metadata=seq.metadata))
    return TabularMSA(seqs)


def hamming_distance_matrix(msa, ignore_sequence_ids=False):
    """Compute Hamming distance matrix of an MSA.

    Parameters
    ----------
    msa : skbio TabularMSA
        Aligned sequences for calculating pairwise Hamming distances
    ignore_sequence_ids : bool
        Default is False. If true, ignore sequence identifier of alignment.
        Useful if identifier got truncated by alignment producing program such
        that different sequences collapse to the same identifier.

    Returns
    -------
    skbio DistanceMatrix
    """
    key = 'id'
    if ignore_sequence_ids:
        key = None
    return DistanceMatrix.from_iterable(msa, hamming, key=key, validate=False)


def cluster_sequences(dm, cutoff):
    """Perform hierarchical clustering based on a distance matrix.

    Parameters
    ----------
    dm : skbio DistanceMatrix
        aligned sequences for calculating pairwise Hamming distances
    cutoff : float
        sequence percent identity cutoff for defining clusters

    Returns
    -------
    list of int
        flat cluster numbers to which sequences are assigned

    Notes
    -----
    returns an empty list if the distance matrix is empty (e.g., there is only
    one input sequence)
    """
    t = 1.0 - cutoff / 100.0
    return list(fcluster(linkage(dm.condensed_form()), t,
                criterion='distance')) if dm.size > 1 else []


def effective_family_size(clusters, length, Nclu=False):
    """Calculate effective family size based on a clustering scheme.

    Parameters
    ----------
    clusters : list of tuple
        clustering scheme (list of (sequence ID, flat cluster number assigned
        to the sequence))
    length : int
        length (number of columns) of the multiple sequence alignment
    Nclu : bool
        return N_cluster instead of N_eff (default: False)

    Returns
    -------
    float or int
        effective family size: N_eff = N_cluster / sqrt(length)
        or number of clusters (N_cluster)
    """
    nc = len(set(clusters))
    return nc if Nclu else nc / math.sqrt(length)


@click.command()
@click.option('--infile', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Input multiple sequence alignment file in A3M format.')
@click.option('--outfile', '-o', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output file of calculated effective family size.')
@click.option('--cutoff', '-c', required=False, type=int, default=100,
              help=('Percent identity cutoff for clustering sequences '
                    '(default: 100).'))
def _calculate_Neff(infile, outfile, cutoff):
    """Parsing arguments for processing.
    """
    msa = parse_msa_file(infile)
    hdm = hamming_distance_matrix(msa)
    clu = cluster_sequences(hdm, cutoff)
    Neff = effective_family_size(clu, msa.shape[1])
    click.echo('Effective family size at %s%% identity: %.3f.'
               % (cutoff, Neff))
    if outfile is not None:
        with open(outfile, 'w') as f:
            f.write('%s\n' % Neff)
    click.echo('Task completed.')


if __name__ == "__main__":
    _calculate_Neff()
