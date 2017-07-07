import os
import numpy as np
import click
from microprot.scripts import process_fasta


@click.command()
@click.option('--infile', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True),
              help='Input protein sequence file in FASTA format.')
@click.option('--outfile', '-o', required=False, default=None,
              type=click.Path(resolve_path=True, readable=True, exists=False),
              help='Output protein sequence file in FASTA format.')
@click.option('--sort_by_len', '-s', required=False, is_flag=True,
              help='Sort sequences by length')
@click.option('--min_len', '-m', required=False, default=1, type=int,
              help='Minimum sequence length to be included in output.')
@click.option('--max_len', '-x', required=False, default=100000, type=int,
              help='Maximum sequence length to be included in output.')
def _process_fasta_input(infile, outfile, sort_by_len, min_len, max_len):
    fp = infile
    fp_name = os.path.splitext(fp)[0]

    fasta = process_fasta.extract_sequences(fp)
    fasta = np.array(fasta)

    suffix = []

    if sort_by_len is True:
        suffix.append('_sorted')
        s = [len(seq) for seq in fasta]
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
        output_fasta = [seq for seq in fasta if min_len <= len(seq) <= max_len]

    # checks if settings changed from default
    if min_len > 1:
        suffix.append('%s%i' % ('_min', min_len))
    if max_len != 100000:
        suffix.append('%s%i' % ('_max', max_len))

    suffix = ''.join(suffix)

    if outfile is None:
        outfile = '%s%s.fasta' % (fp_name, suffix)
    process_fasta.write_sequences(output_fasta, outfile)


if __name__ == "__main__":
    _process_fasta_input()
