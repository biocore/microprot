#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Tomasz Kosciolek"
__version__ = "1.01b"
__last_update__ = "12/05/2016"


import click
import sys
from skbio import Sequence


def extract_sequences(infile, outfile, seqidx=0):
    """ Extract sequence(s) from a multi-sequence FASTA file

    Parameters
    ----------
    infile : str
        file path to input multi-sequence FASTA file
    outfile : str
        file path to output FASTA file
    seqidx : int (optional)
        n-th sequence of the input file to extract
        (default: 0, standing for all sequences)
    """


@click.command()
@click.option('--infile', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
                              file_okay=True),
              help='Input protein sequence file in FASTA format')
@click.option('--outfile', '-o', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=False,
                              file_okay=True),
              help='Output protein sequence file in FASTA format')
@click.option('--seqidx', '-s', default=0,
              help='Extract n-th sequence (default: all sequences)')
def _main(infile, outfile, seqidx):
    """ Parsing arguments for pre-processing
    """
    extract_sequences(infile, outfile, seqidx)


if __name__ == "__main__":
    main()
