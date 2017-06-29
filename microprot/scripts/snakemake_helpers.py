import os
from datetime import datetime
from microprot.scripts import process_fasta
from skbio import io


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
    if msa_ext != '.a3m':
        msa_ext = '.a3m'
    with open(''.join([msa_dir, msa_ext]), 'r') as f:
        lines = f.readlines()
    msa_size = len([line for line in lines if line.startswith('>')])-1
    return msa_size


def parse_inputs(inp_file=None, inp_path=None, inp_from=None, inp_to=None,
                 microprot_inp=None, microprot_out=None):
    SEQ_file = inp_file
    SEQ_path = inp_path
    SEQ_from = inp_from
    SEQ_to = inp_to
    SEQS = process_fasta.extract_sequences('%s/%s' % (SEQ_path, SEQ_file),
                                           identifiers=(SEQ_from, SEQ_to))
    SEQ_ids = []
    for SEQ in SEQS:
        _seq = SEQ.metadata['id']
        _seq = _seq.replace('/', '_')
        _seq = _seq.replace('\\', '_')
        _seq = _seq.replace('|', '_')
        SEQ_ids.append(_seq)
        SEQ.metadata['id'] = _seq
        io.write(SEQ, format='fasta', into='%s/%s.fasta' % (microprot_inp,
                                                            _seq))
        io.write(SEQ, format='fasta',
                 into='%s/%s' % (microprot_out, 'processed_sequences.fasta'))
        return SEQ_ids

# TODO
# adding MSA size information omitted for now
def append_db(fname, step=None, version=1, db_fp='/tmp/protein_db.index'):
    prots = process_fasta.extract_sequences(fname)
    # fp = trim(fname, '/')
    for prot in prots:
        prot_name = prot.metadata['id']
        timestamp = str(datetime.now())
        """
        `msa_fp` needs to be retrieved from 1 level up (e.g. 01 folder,
        instead of 02 folder. For future consideration. Added as issue #48
        """
        # msa_fp = '%s/%s' % (fp, prot_name)
        # msa_size = msa_size(msa_fp)

        # > protein_name # source # commit_no # timestamp # (msa_size)
        append = '> %s # %s # %i # %s\n' % (prot_name,
                                            step,
                                            version,
                                            timestamp)
        with open(db_fp, 'a') as f:
            f.write(append)
