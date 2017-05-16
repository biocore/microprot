import os
import urllib
import click
from xml.dom import minidom
from skbio import io

"""
Notes:
This script has two functions:
1. Extract select protein(s) from a multi-sequence Fasta file:
    python processing.py -i input.faa -d 10 -o output.faa
    python processing.py -i input.faa -d 1K5N.B -o output.faa
    python processing.py -i input.faa -d 1,2,3 -o output.faa
    python processing.py -i input.faa -d list.txt -o output.faa
2. Extract representative proteins (local or remote) from PDB:
    python processing.py -i pdb_seqres.txt -r represents.xml -o output.faa
Or:
    python processing.py -i pdb_seqres.txt -r 100 output.faa
Stats:
    pdb_seqres.txt has 381126 sequences
    representatives?cluster=100 has 70132 protein IDs, as of Jan. 16, 2017.
"""


def extract_sequences(infile, identifiers=None):
    """ Extract sequence(s) from a multi-sequence FASTA file.

    Parameters
    ----------
    infile : str
        file path to input multi-sequence FASTA file
    identifiers : str or list of str (optional)
        IDs (names) or indexes (n-th sequence in the file) of the proteins to
        extract, may be:
        1) Python list
        2) comma-separated list
        3) file path to the list (one ID or index per line)
        if omitted, all sequences will be extracted

    Return
    ------
    list of skbio Sequence
        extracted protein sequences
    """
    l, ids, indexes = [], set(), set()
    if identifiers:
        if type(identifiers) is list:  # Python list
            l = identifiers
        elif os.path.isfile(identifiers):  # read from a local file
            with open(identifiers, 'r') as f:
                l = f.read().splitlines()
        else:  # comma-separated list
            l = list(map(str.strip, identifiers.split(',')))
        for i in l:
            if i.isdigit():
                indexes.add(int(i))  # index of this protein in the file
            else:
                ids.add(i)  # protein ID (name)
    seqs = []
    for i, seq in enumerate(io.read(infile, format='fasta')):
        if ids:
            if seq.metadata['id'] in ids:
                seqs.append(seq)
        elif indexes:
            if i+1 in indexes:  # indexes start with 1, not 0
                seqs.append(seq)
        else:
            seqs.append(seq)
    return seqs


def write_sequences(seqs, outfile):
    """ Write sequence(s) into a multi-sequence FASTA file.

    Parameters
    ----------
    seqs : str
        list of skbio Sequence
    outfile : str
        file path to output multi-sequence FASTA file
    """
    def outseqs():
        for seq in seqs:
            yield seq
    io.write(outseqs(), format='fasta', into=outfile)


def read_representatives(represent):
    """ Read representative protein IDs from file or server

    Parameters
    ----------
    represent : str
        file path to list of representative protein IDs in XML format, or an
        integer representing the clustering level by which the list will be
        downloaded from the PDB server

    Return
    ------
    list of str
        protein IDs

    Raises
    ------
    ValueError
        if neither repfile nor cluster is specified
    """
    ids = set()
    xmldoc = None
    if os.path.isfile(represent):  # read from a local file
        xmldoc = minidom.parse(represent)
    elif represent.isdigit():  # download from the PDB server
        url_str = ('http://www.rcsb.org/pdb/rest/representatives?cluster=%s'
                   % represent)
        xml_str = urllib.request.urlopen(url_str).read()
        xmldoc = minidom.parseString(xml_str)
    else:
        raise ValueError('Error: You must specify a local file path or the '
                         'clustering level.')
    for x in xmldoc.getElementsByTagName('pdbChain'):
        l = x.attributes['name'].value.split('.')
        ids.add('%s_%s' % (l[0].lower(), l[1]))
    return sorted(list(ids))


@click.command()
@click.option('--infile', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Input protein sequence file in FASTA format.')
@click.option('--outfile', '-o', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output protein sequence file in FASTA format.')
@click.option('--identifiers', '-d', required=False,
              help='Indexes or IDs of protein sequences to extract. May be a '
              'comma-separated list, or an external file with one entry per '
              'line. If not specified, all sequences will be extracted.')
@click.option('--represent', '-r', required=False,
              help='Generate representative protein sequences. One may either '
              'provide a local protein ID list file (XML format), or let the '
              'program download it from the PDB server by specifying a '
              'desired clustering level.')
def _processing(infile, outfile, identifiers, represent):
    """ Parsing arguments for processing
    """
    if represent:
        identifiers = read_representatives(represent)
        click.echo('Number of representative proteins: %s' % len(identifiers))
    seqs = extract_sequences(infile, identifiers)
    click.echo('Number of extracted proteins: %s' % len(seqs))
    if seqs:
        write_sequences(seqs, outfile)
    click.echo('Task completed.')


if __name__ == "__main__":
    _processing()
