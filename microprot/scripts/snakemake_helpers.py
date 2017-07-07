import os
import textwrap
from datetime import datetime
from microprot.scripts import process_fasta
from skbio import io


def not_empty(fname):
    """ Check if a directory exists and/or file exists and is not empty
    Parameters
    ----------
    fname: string
        file path with filename

    Returns
    -------
    bool
    """
    if os.path.exists(fname) and os.path.getsize(fname) > 0:
        return True
    else:
        return False


def trim(inp_str, symbol):
    """ Trim a string up to a symbol
    Parameters
    ----------
    inp_str : string
        input string
    symbol : string
        symbol to split on

    Returns
    -------
    out : string
        trimmed output
    """
    if isinstance(inp_str, str) and isinstance(symbol, str):
        out = symbol.join(inp_str.split(symbol)[:-1])
    else:
        raise TypeError('Trim function requires strings as input!')
    return out


def msa_size(msa_fp):
    msa_dir, msa_ext = os.path.splitext(os.path.abspath(msa_fp))
    if msa_ext != '.a3m':
        msa_ext = '.a3m'
    with open(''.join([msa_dir, msa_ext]), 'r') as f:
        lines = f.readlines()
    msa_size = len([line for line in lines if line.startswith('>')])-1
    return msa_size


def parse_inputs(inp_fp=None, inp_from=None, inp_to=None,
                 microprot_inp=None, microprot_out=None):

    for _dir in [microprot_inp, microprot_out]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    SEQS = process_fasta.extract_sequences(inp_fp,
                                           identifiers=(inp_from, inp_to))
    SEQ_ids = []
    processed_fh = open('%s/%s' % (microprot_out,
                                   'processed_sequences.fasta'), 'a')
    for i, SEQ in enumerate(SEQS):
        _seq = SEQ.metadata['id']
        _seq = _seq.replace('/', '_')
        _seq = _seq.replace('\\', '_')
        _seq = _seq.replace('|', '_')
        SEQ_ids.append(_seq)
        SEQ.metadata['id'] = _seq
        io.write(SEQ, format='fasta', into='%s/%s.fasta' % (microprot_inp,
                                                            _seq))
        io.write(SEQ, format='fasta',
                 into=processed_fh)
    processed_fh.close()
    return SEQ_ids


# TODO
# adding MSA size information omitted for now
def write_db(fname, step=None, version=1, db_fp='/tmp/protein_db'):
    prots = process_fasta.extract_sequences(fname)
    # fp = trim(fname, '/')
    for prot in prots:
        prot_name = prot.metadata['id']
        timestamp = str(datetime.now()).split('.')[0]
        """
        `msa_fp` needs to be retrieved from 1 level up (e.g. 01 folder,
        instead of 02 folder. For future consideration. Added as issue #48
        """
        # msa_fp = '%s/%s' % (fp, prot_name)
        # msa_size = msa_size(msa_fp)

        # > protein_name # source # commit_no # timestamp # (msa_size)
        append_idx = '>%s # %s # %i # %s\n' % (prot_name,
                                               step,
                                               version,
                                               timestamp)
        with open('%s.index' % db_fp, 'a') as f:
            f.write(append_idx)

        append_seq = '>%s\n%s\n' % (prot_name, textwrap.fill(str(prot[:]), 70))
        with open(db_fp, 'a') as f:
            f.write(append_seq)


def append_db(source_root, dest_root):
    for ext in ['', '.index']:
        _dest = open('%s%s' % (dest_root, ext), 'a')
        _source_path = '%s%s' % (source_root, ext)
        if not_empty(_source_path):
            _dest.write(open(_source_path, 'r').read())
