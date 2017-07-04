import numpy as np
import click
import os
from skbio import io
from microprot.scripts import process_fasta


@click.command()
@click.option('--infile', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Input protein sequence file in FASTA format.')
@click.option('--outfile', '-o', required=False,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output protein sequence file in FASTA format.')
@click.option('--sort', '-s', required=False, is_flag=True,
              help='Sort sequences by length')
@click.option('--min_len', '-m', required=False, default=1,
              help='Minimum sequence length to be included in output.')
@click.option('--max_len', '-x', required=False, default=10000,
              help='Maximum sequence length to be included in output.')
def _process_fasta_input(infile, outfile, sort, min_len, max_len):
    fp = infile
    fp_name = os.path.splitext(fp)[0]

    fasta = process_fasta.extract_sequences(fp)
    fasta = np.array(fasta)

    suffix = []

    if sort is True:
        suffix.append('_sorted')
        s = []
        for seq in fasta:
            s.append(len(seq))
        idx = sorted(range(len(s)), key=lambda k: s[k])
        fasta = fasta[idx]
        for i, seq in enumerate(fasta):
            if len(seq) >= min_len:
                output_fasta = fasta[i:]
                break
        for i, seq in enumerate(output_fasta):
            if len(seq) > max_len:
                output_fasta = output_fasta[:i]
                break
    else:
        idx = []
        for i, seq in enumerate(fasta):
            if len(seq) >= min_len and len(seq) <= max_len:
                idx.append(i)
        output_fasta = fasta[idx]

    if min_len > 1:
        suffix.append('%s%i' % ('_min', min_len))
    if max_len != 10000:
        suffix.append('%s%i' % ('_max', max_len))

    suffix = ''.join(suffix)

    if outfile:
        process_fasta.write_sequences(output_fasta, outfile)
    else:
        process_fasta.write_sequences(output_fasta, '%s%s.fasta' % (fp_name,
                                                                    suffix))

if __name__ == "__main__":
    _process_fasta_input()
