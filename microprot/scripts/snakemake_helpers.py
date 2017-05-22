import os
from datetime import datetime
from microprot.scripts import process_fasta


def not_empty(fname):
    if os.path.exists(fname) and os.path.getsize(fname) > 0:
        return True
    else:
        return False


def trim(inp_str, symbol):
    if isinstance(inp_str, str) and isinstance(symbol, str):
        out = symbol.join(inp_str.split(symbol)[:-1])
    else:
        raise TypeError('Trim function requires strings as input!')
    return out


def msa_size(msa_fp):
    msa_dir, msa_ext = os.path.splitext(os.path.abspath(msa_fp))
    if msa_ext is not '.a3m':
        msa_ext = '.a3m'
    with(''.join([msa_dir, msa_ext]), 'r') as f:
        lines = f.readlines()
    msa_size = len([line for line in lines if line.startswith('>')])-1
    return msa_size


def append_db(fname, step=None, version=1, db_fp='/tmp/protein_db.index'):
    prots = process_fasta.extract_sequences(fname)
    for prot in prots:
        prot_name = prot.metadata['id']
        timestamp = str(datetime.now())
        msa_size = msa_size(fname)

        # > protein_name # source # commit_no # timestamp # MSA_size
        append = '> %s # %s # %i # %s # %i\n' % (prot_name,
                                                 step,
                                                 version,
                                                 timestamp,
                                                 msa_size)
        with open(db_fp, 'a') as f:
            f.write(append)


def search_X(seq, split, output_base, shell_func):
    for d_boundaries, domain in split_sequence(seq, split):
        n = '%s-%s' % (str(d_boundaries[0]), str(d_boundaries[1]))
        out = '%s_%s.%s' % (output_base, n, 'out')
        a3m = '%s_%s.%s' % (output_base, n, 'a3m')
        with open(output_base, 'a') as f:
            f.write(out)
            f.write(a3m)
        shell_func
        pass
