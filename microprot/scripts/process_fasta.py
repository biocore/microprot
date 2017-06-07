import os
import urllib
import click
from xml.dom import minidom
from skbio import io
from skbio.sequence import Sequence


"""
Notes:
This script has three functionalities:
1. Extract select protein(s) from a multi-sequence Fasta file:
    python process_fasta.py -i input.faa -d 10 -o output.faa
    python process_fasta.py -i input.faa -d 1K5N.B -o output.faa
    python process_fasta.py -i input.faa -d 1,2,3 -o output.faa
    python process_fasta.py -i input.faa -d 1..3 -o output.faa
    python process_fasta.py -i input.faa -d list.txt -o output.faa
2. Extract representative proteins (local or remote) from PDB:
    python process_fasta.py -i pdb_seqres.txt -r represents.xml -o output.faa
Or:
    python process_fasta.py -i pdb_seqres.txt -r 100 output.faa
3. python process_fasta.py -i input.faa --split -o newdir
Stats:
    pdb_seqres.txt has 381126 sequences
    representatives?cluster=100 has 70132 protein IDs, as of Jan. 16, 2017.
"""


def extract_sequences(infile, identifiers=None):
    """Extract sequence(s) from a multi-sequence FASTA file.

    Parameters
    ----------
    infile : str
        file path to input multi-sequence FASTA file
    identifiers :
        int
            sequence index (n-th sequence in the file)
        str
            sequence ID (name) or index
                numeric str is treated as index instead of ID
            comma-separated sequence IDs or indexes
            file path to sequence list (one ID or index per line)
            sequence index range as "start..end" (both included)
                start must be smaller or equal to end
        list of int
            sequence indexes
        list of str
            sequence IDs or indexes
        tuple of two int's
            sequence index range as (start, end)
        if omitted, all sequences will be extracted

    Returns
    -------
    list of skbio Sequence
        extracted protein sequences

    Raises
    ------
    ValueError
        if tuple (index range) is not in (start, end) form
        if index range str is not formatted as "start, end"
        if the data type of identifiers is incorrect
    """
    l, ids, indexes = [], set(), set()
    if identifiers:
        # IDs or indexes as list
        if isinstance(identifiers, list):
            l = identifiers
        # start and end indexes as tuple of int
        elif isinstance(identifiers, tuple):
            if len(identifiers) == 2 \
                and all(isinstance(n, int) for n in identifiers) \
                    and 0 < identifiers[0] <= identifiers[1]:
                l = list(range(identifiers[0], identifiers[1] + 1))
            else:
                raise ValueError('Error: Index range must be a tuple of '
                                 '(start, end).')
        elif isinstance(identifiers, str):
            # read from a file
            if os.path.isfile(identifiers):
                with open(identifiers, 'r') as f:
                    l = f.read().splitlines()
            # start and end indexes as str
            elif '..' in identifiers:
                l = identifiers.split('..')
                if len(l) == 2 \
                    and all(n.isdigit() for n in l) \
                        and 0 < int(l[0]) <= int(l[1]):
                    l = list(range(int(l[0]), int(l[1]) + 1))
                else:
                    raise ValueError('Error: Index range must be formatted as '
                                     '"start..end".')
            # IDs or indexes as str (single or comma-separated list)
            else:
                l = list(map(str.strip, identifiers.split(',')))
        # index as int
        elif isinstance(identifiers, int):
            l = [identifiers]
        else:
            raise ValueError('Error: Incorrect data type of identifiers.')
        for i in l:
            if isinstance(i, int):  # index of this protein in the file
                indexes.add(i)
            elif i.isdigit():
                indexes.add(int(i))
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
    """Write sequence(s) into a multi-sequence FASTA file.

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
    """Read representative protein IDs from file or server.

    Parameters
    ----------
    represent : str
        file path to list of representative protein IDs in XML format, or an
        integer representing the clustering level by which the list will be
        downloaded from the PDB server

    Returns
    -------
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


def split_fasta(seqs, prefix=None, outdir=None):
    """Split a multi-protein FASTA file into single-sequence FASTAs.

    Parameters
    ----------
    seqs : list of skbio.sequence
        List of skbio protein sequences with a protein name in header
    prefix : str
        Name prefix to be added to output single-sequence FASTA files

    Raises
    ------
    TypeError
        seqs needs to be a filepath or skbio.sequence object
    """

    if isinstance(seqs, str):
        if os.path.exists(seqs):
            seqs = extract_sequences(seqs)
        else:
            raise TypeError('split_fasta sequence input is not a filepath or '
                            'filepath does not exist.')
    elif isinstance(seqs, list):
        if len(seqs) == 0:
            raise ValueError('Empty list provided to split_fasta')
        else:
            if not isinstance(seqs[0], Sequence):
                raise TypeError('Object you provided to split_fasta is not '
                                'a skbio.sequence object')
    else:
        raise TypeError('split_fasta input sequences need to be a filepath or'
                        'skbio.sequence object.')

    if not outdir:
        outdir = os.getcwd()
    elif not os.path.exists(outdir):
        os.makedirs(outdir)
    for seq in seqs:
        if prefix:
            io.write(seq, format='fasta', into='%s/%s_%s.fasta' %
                     (outdir, prefix, seq.metadata['id']))
        else:
            io.write(seq, format='fasta', into='%s/%s.fasta' %
                     (outdir, seq.metadata['id']))


@click.command()
@click.option('--infile', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Input protein sequence file in FASTA format.')
@click.option('--outfile', '-o', required=False,
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
@click.option('--split', '-s', required=False, is_flag=True,
              help='Split multi-sequence FASTA file into single sequence FAST'
                   'As with sequence name as file name.')
@click.option('--prefix', '-p', required=False, default=None,
              help='Prefix to be used in splitting multi-sequence FASTA file.')
def _process_fasta(infile, outfile, identifiers, represent, split, prefix):
    """Parsing arguments for processing.
    """
    if not outfile and not split:
        raise IOError('No outfile or split flag used. You need to specify at '
                      'least one.')

    if represent:
        identifiers = read_representatives(represent)
        click.echo('Number of representative proteins: %s' % len(identifiers))
    seqs = extract_sequences(infile, identifiers)
    click.echo('Number of extracted proteins: %s' % len(seqs))
    if split:
        split_fasta(seqs, prefix, outfile)
    elif seqs:
        write_sequences(seqs, outfile)
    click.echo('Task completed.')


if __name__ == "__main__":
    _process_fasta()
