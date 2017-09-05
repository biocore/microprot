import sys
import re

import click
from skbio import Protein

_HEADER = ['No',
           'Hit',
           'Prob',
           'E-value',
           'P-value',
           'Score',
           'SS',
           'Cols',
           'Query HMM',
           'Template HMM',
           '#match states']


def is_overlapping(intA, intB):
    ''' Test if two intervals overlap.

    Parameters
    ----------
    intA : (int, int)
        A pair of coordinates, e.g. (20, 100). Left component must be <= right.
    intB : (int, int)
        A pair of coordinates, e.g. (50, 150). Left component must be <= right.

    Returns
    -------
    Bool. True if intA overlaps with intB. False otherwise.

    Raises
    ------
    ValueError
        If left component of intA or intB is greater than its right component.
    '''
    for x in zip(["intA", "intB"], [intA, intB]):
        if (x[1][0] > x[1][1]):
            raise ValueError("Left component of %s cannot be larger than its "
                             "right component, i.e. the interval must be "
                             "forward oriented." % x[0])

    return not (((intA[0] > intB[0]) & (intA[0] > intB[1])) |
                ((intA[1] < intB[0]) & (intA[1] < intB[1])))


def _parse_hit_summary_line(line):
    """ Parses a single hit summary line of HHsearch.

    Fields are described in [1]. Unfortunately, there is no single dedicated
    separator char plus the second field might contain whitespaces. Thus, I
    split the line on \s and consume the first element for the 'hit index'. The
    'description' are all elements that are not consumed by the last 9 fields.

    Parameters
    ----------
    line: str
        The line that should be parsed and contain the expected information.

    Returns
    -------
    A dict holding the information of the HHsearch line for a hit.

    Raises
    ------
    ValueError
        If some fields cannot be casted to int or float, which might be an
        indication that the provided line is no HHsearch hit summary line.

    [1] Johannes Soeding, "Quick guide to HHsearch"
        ftp://ftp.tuebingen.mpg.de/pub/protevo/HHsearch/HHsearch1.5.01/HHsearc
        h-guide.pdf
    """

    fields = line.rstrip().split()

    # In some cases the field for Template HMM start-end positions is too long,
    # such that no space is left for the delimiting char to the neighboring
    # column (length). We detect those pathological cases by testing if the
    # last field does not start with '('. We then need to 'manually' split this
    # field. See issue https://github.com/biocore/microprot/issues/33
    if not fields[-1].startswith('('):
        pos, length = fields[-1].split('(')
        del fields[-1]
        fields.extend([pos, '(' + length])

    hit = {}

    try:
        # Index of hit
        hit[_HEADER[0]] = int(fields[0])

        # First 30 characters of domain description (from nameline of query
        # sequence)
        hit[_HEADER[1]] = " ".join(fields[1:-9])

        # Probability of target to be a true positive For the probability of
        # being a true positive, the secondary structure score in column 7 is
        # taken into account, together with the raw score in column 6
        # (’Score’). True positives are defined to be either globally
        # homologous or they are at least locally similar in structure. More
        # precisely, the latter criterion demands that the MAXSUB score between
        # query and hit is at least 0.1. In almost all cases the structural
        # similarity will we be due to a global OR LOCAL homology between query
        # and target.
        hit[_HEADER[-9]] = float(fields[-9])

        # Expect-value E-value and P-value are calculated without taking the
        # secondary structure into account! The E-value gives the average
        # number of false positives (’wrong hits’) with a score better than the
        # one for the target when scanning the datbase. It is a measure of
        # reliability: E-values near to 0 signify a very reliable hit, an
        # E-value of 10 means about 10 wrong hits are expected to be found in
        # the database with a score at least this good.
        hit[_HEADER[-8]] = float(fields[-8])

        # P-value The P-value is just the E-value divided by the number of
        # sequences in the database. It is the probability that in a PAIRWISE
        # comparison a wrong hit will score at least this good.
        hit[_HEADER[-7]] = float(fields[-7])

        # Raw score, does not include the secondary structure score
        hit[_HEADER[-6]] = float(fields[-6])

        # Secondary structure score This score tells how well the PSIPRED-
        # predicted (3-state) or actual DSSP-determined (8-state) secondary
        # structure sequences agree with each other. PSIPRED confidence values
        # are used in the scoring, low confidences getting less statistical
        # weight.
        hit[_HEADER[-5]] = float(fields[-5])

        # The number of aligned Match columns in the HMM-HMM alignment.
        hit[_HEADER[-4]] = int(fields[-4])

        # Range of aligned match states from query HMM
        hit[_HEADER[-3]] = fields[-3]

        # Range of aligned match states from target HMM
        hit[_HEADER[-2]] = fields[-2]

        # Number of match states in target HMM
        hit[_HEADER[-1]] = int(fields[-1][1:-1])

        return hit
    except ValueError:
        raise ValueError("Unexpected field. Check if line is a HHsearch hit"
                         " summary line.")


def _parse_hit_block(block):
    """ Parse one alignment block of HHsearch output.

    Unfortunately, there is no strict format description. It is unclear which
    lines might appear. Identifier might have white spaces. I assume that there
    will always be a line starting with 'Q Consensus', having some integer to
    describe the start position of the query, followed by a part of the aligned
    query sequence. Delimiter are unclear, I guess it is \s+.
    Due to all those uncertainties the code might look confusing at first
    sight.

    Parameters
    ----------
    block : str
        The lines of the HHsearch output regarding to one alignment of query
        sequence and target HMM. It might be wrapped into several parts.

    Returns
    -------
    A dict holding all information about the aligment.
    """

    hit = {}
    lines = block.split('\n')

    # Index of hit
    hit['No'] = int(lines[0].split()[1])

    # description of domain
    hit['Hit'] = lines[1][1:]

    # line with stats
    for field in lines[2].rstrip().split('  '):
        key, value = field.split('=')
        hit[key] = value
    hit['Similarity'] = float(hit['Similarity'])
    hit['Score'] = float(hit['Score'])
    hit['Probab'] = float(hit['Probab'])
    hit['Identities'] = float(hit['Identities'][:-1])/100
    hit['Aligned_cols'] = int(hit['Aligned_cols'])
    hit['E-value'] = float(hit['E-value'])
    hit['Sum_probs'] = float(hit['Sum_probs'])
    hit['Template_Neff'] = float(hit['Template_Neff'])

    idx = 4
    block = {}
    while(idx+1 < len(lines)):
        # determin start and end column of alignment content
        [coord, content] = lines[idx+1].split()[2:4]
        coord_pos = lines[idx+1].find(coord)
        startCol = lines[idx+1].find(content, coord_pos)
        endCol = startCol + len(content)
        while((idx+1 < len(lines)) & (lines[idx] != '')):
            start = None
            if lines[idx].startswith(' '):
                name = 'column score'
            else:
                parts = lines[idx][:startCol].split()
                try:
                    start = int(parts[-1])  # check if last part is a number
                    name = lines[idx][:lines[idx][:startCol].find(parts[-1])]
                except ValueError:
                    name = lines[idx][:startCol]
                name = name.rstrip()
            if name not in block:
                block[name] = {'sequence': lines[idx][startCol:endCol]}
                if start is not None:
                    block[name]['start'] = start
            else:
                block[name]['sequence'] += lines[idx][startCol:endCol]
            if lines[idx][endCol:] != '':
                fields = lines[idx][endCol:].split()
                block[name]['end'] = int(fields[0])
                block[name]['totallen'] = int(fields[1][1:-1])
            idx += 1
        idx += 2

    hit['alignment'] = block

    return hit


def parse_pdb_match(filename):
    """ Parse an HHsearch output file.

    Parameters
    ----------
    filename : str
        Path to the HHsearch output file that should be parsed.

    Returns
    -------
    A list of HHsearch hits. Each hit is a dict, holding all its information.

    Raises
    ------
    IOError
        If the file cannot be read.
    """
    hits = []
    try:
        fh = open(filename, 'r')
        line = ""

        # read until header of summary table is found
        while(len(set(line.rstrip().split()) & set(_HEADER)) < 8):
            line = fh.readline()
        line = fh.readline()  # summary table header line
        # read all summary lines
        while(line != '\n'):
            hits.append(_parse_hit_summary_line(line))
            line = fh.readline()

        # read the alignments
        reachedEOF = False
        alignments = []
        block = ''
        while(not reachedEOF):
            line = fh.readline()
            reachedEOF = line == ''
            if line.startswith('No ') & (len(line.split()) == 2):
                if block != '':
                    alignments.append(_parse_hit_block(block))
                block = line
            else:
                block += line
        if block != '':
            alignments.append(_parse_hit_block(block))

        # merge hit summary and according alignment
        for idx in range(len(hits)):
            for key in alignments[idx]:
                hits[idx][key] = alignments[idx][key]
            # duplicate information, but more precicse on other place
            del hits[idx][_HEADER[2]]  # Prob
            del hits[idx][_HEADER[8]]  # Query HMM
            del hits[idx][_HEADER[9]]  # Template HMM
            del hits[idx][_HEADER[10]]  # match states

        fh.close()
        return hits
    except IOError:
        raise IOError('Cannot read file "%s"' % filename)


def select_hits(hits, e_value_threshold=0.001):
    """ Picking HHsearch hits from the list of all hits.

    We take those hits that a) have an e-Value <= e_value_threshold and b)
    does not overlap with any previously picked hit.

    Parameters
    ----------
    hits: [dict]
        The sorted list of hits from an HHsearch.
    e_value_threshold: float
        The threshold for a maximal e-Value for a hit to pick. Default: 0.001

    Returns
    -------
    [dict] a potentially smaller list of HHsearch hits.
    """
    good_hits = []

    for hit in hits:
        if hit['E-value'] < e_value_threshold:
            id_hit = get_q_id(hit)
            # check if there is an overlap to previous results
            noOverlap = True
            for good in good_hits:
                id_good = get_q_id(good)
                if is_overlapping((hit['alignment'][id_hit]['start'],
                                   hit['alignment'][id_hit]['end']),
                                  (good['alignment'][id_good]['start'],
                                   good['alignment'][id_good]['end'])):
                    noOverlap = False
                    break
            if noOverlap:
                good_hits.append(hit)

    return good_hits


def report_hits(hits):
    """ Return condensed information about HHsearch hits.

    Parameters
    ----------
    hits : [dict]
        A list of HHsearch hits.

    Returns
    -------
    For each hit: a) the pdb_id, b) the protein sub-sequence of the query
    covered by this hit, c) the start position of the hit relative to the query
    and d) the end position of the hit relative to the query (first AA in query
    is position 1, not 0).
    """
    report = []
    for hit in hits:
        # find the right key for the query sequence information
        q_keys = get_q_id(hit)

        # compose a new dict as report for this hit.
        report.append({'pdb_id': hit['Hit'].split()[0],
                       'covered_sequence': hit['alignment'][q_keys]
                       ['sequence'].replace('-', ''),
                       'start': hit['alignment'][q_keys]['start'],
                       'end': hit['alignment'][q_keys]['end']})
    return report


def report_uncovered_subsequences(hits, query, min_subseq_len=40):
    """ Returns sub-sequences of query that is not covered by hits.

    Parameters
    ----------
    hits : [dict]
        A list of HHsearch hits.
    query : str
        The protein sequence that was used to generate the HHsearch hits, i.e.
        the content of the input file for HHsearch.
    min_subseq_len : int
        A threhold for the minimal sub-sequence length that should be reported.
        Default = 40

    Returns
    -------
    A list of dicts, where each entry corresponds to a sub-sequence of the
    query that is not covered by the HHsearch hits and whose length is at least
    min_subseq_len. The dict holds the keys a) 'sequence' for the plain AA
    sequence, b) 'start' for the start position of the reported sub-sequence,
    relative to the query sequence and c) 'end' the relative end position.
    (positions start with 1, not with 0)
    """
    covered = []
    for hit in hits:
        # find the right key for the query sequence information
        q_keys = get_q_id(hit)
        covered.append((hit['alignment'][q_keys]['start'],
                        hit['alignment'][q_keys]['end']))

    uncovered = []
    curEnd = 0
    for sub in sorted(covered):
        unc = {'sequence': query[curEnd:sub[0]-1],
               'start': curEnd+1,
               'end': sub[0]-1}
        uncovered.append(unc)
        curEnd = sub[1]
    unc = {'sequence': query[curEnd:],
           'start': curEnd+1,
           'end': len(query)}
    uncovered.append(unc)

    return [u for u in uncovered if (u['sequence'] != '') &
                                    (len(u['sequence']) >= min_subseq_len)]


def get_q_id(hit):
    """ Returns the query ID for a hit.

    Parameters
    ----------
    A hit parsed from an HHsearch output file, i.e. dict with at least the key
    'alignment' which is a dict by itself and comes at least with the key
    'Q xxx' where xxx is some identifier. The value for this 'Q xxx' key is a
    third dict which needs to contain the key 'sequence'.

    Returns
    -------
    str : The query ID starting with 'Q '.

    Notes
    -----
    Each 'hit' has one 'alignment', which comes with different lines. One of
    those lines is 'Q consensus'. Another line is called 'Q xxx' where xxx is
    the ID of the input query sequence. This function find this 'Q xxx' ID.
    We assume that there are only two line names starting with 'Q '.
    """
    # find the right ID
    _id = [_id for _id in hit['alignment'].keys()
           if _id.startswith('Q') and
           _id != 'Q Consensus'][0]
    return _id


def frag_size(hit):
    """ Compute the fragment length of a hit.

    Parameters
    ----------
    A hit parsed from an HHsearch output file, i.e. dict with at least the key
    'alignment' which is a dict by itself and comes at least with the key
    'Q xxx' where xxx is some identifier. The value for this 'Q xxx' key is a
    third dict which needs to contain the key 'sequence'.

    Returns
    -------
    The length of the un-gapped sequence for the given hit."""
    # find the right ID
    _id = get_q_id(hit)
    subseq = hit['alignment'][_id]['sequence']
    # remove gap characters and return length
    return len(subseq) - subseq.count('-')


def correct_header_positions(header, new_start, new_end, out=sys.stderr):
    """Prevents chaining of start-end positions.

    Parameters
    ----------
    header : str
        The existing fasta header string with start-end positions as a suffix.
        Valid header format is specified as the regular expression:
        '(.+)_(\d+)\-(\d+)$'
        Examples for valid header: 'test_1_251-330', 'NZ_GG666849.1_2_251-330'
        Examples for invalid headers: 'NZ_GG666849.1_2', 'test_1_251:330'
    new_start : int
        The new start position, relative to the interval found in the header.
        Example: header 'test_1_251-330', new_start 50 -> updated_start 300
    new_end : int
        The new end position, relative to the interval found in the header.
        Example: header 'test_1_251-330', new_end 70 -> updated_start 320
    out : file handle
        File handle into which verbosity information should be printed.

    Returns
    -------
    str : the new fasta header with updated start, end positions.

    Raises
    ------
    ValueError if
        a) input header string does not follow expected format
        b) old and new start, end positions are incompatible.
    """
    pattern = re.compile('(.+)_(\d+)\-(\d+)$')

    hit = pattern.match(header)
    if hit is None:
        return '%s_%i-%i' % (header, new_start, new_end)
    else:
        prefix, old_start, old_end = hit.group(1), hit.group(2), hit.group(3)
        len_new = int(new_end) - int(new_start) + 1
        len_old = int(old_end) - int(old_start) + 1
        if len_new > len_old:
            raise ValueError("Intervals are not nested!")
        if new_start > len_old:
            raise ValueError("New interval out of old interval borders!")

        updated_start = int(old_start) + int(new_start) - 1
        updated_end = int(old_start) + int(new_end) - 1

        return "%s_%i-%i" % (prefix, updated_start, updated_end)


def mask_sequence(hhsuite_fp, fullsequence_fp, subsequences_fp=None,
                  min_prob=None, max_pvalue=None, max_evalue=None,
                  min_fragment_length=0):
    """ Splits a protein sequence according to HHsuits results.

    The returned sub-sequences will seamlessly build the full sequence if
    re-concatenated.

    Parameters
    ----------
    hhsuite_fp : str
        Filepath to HHblits/HHsearch output.
    fullsequence_fp : str
        Filepath to the protein sequence of the original query.
    subsequences_fp : str
        Filepath to which sub-sequences are written as a multiple fasta file.
        Each sequence makes up one header and one sequence file, i.e. sequences
        are not wrapped.
        Two files will be produced, suffixed by '.match' and '.non_match'. The
        first holds sub-sequences of hits, the second holds the none-hit
        covered subsequences.
        Default: None, i.e. no file is written.
    min_prob: float
        Minimal probability of a hit to be included in the resulting list.
        Note: probabilities are in the range of 100.0 to 0.0.
        Default: None, i.e. no filtering on probability.
    max_pvalue: float
        Maximal P-value of a hit to be included in the resulting list.
        Default: None, i.e. no filtering on P-value.
    max_evalue: float
        Maximal E-value of a hit to be included in the resulting list.
        Default: None, i.e. no filtering on E-value.
    min_fragment_length: int
        Minimal fragment length of a hit to be included in the resulting list.
        Default: 0, i.e. no filtering on fragment length.

    Returns
    -------
    [(str, str)] where first component is a fasta header, while the second is
    its fasta sequence.

    Raises
    ------
    IOError
        If the file cannot be written.

    Notes
    -----
    A hit must satisfy ALL filtering options (min_prob, max_pvalue, max_evalue,
    min_fragment_length) to be included in the resulting list.
    """

    # parse hits from file
    hits = parse_pdb_match(hhsuite_fp)

    # filter hits
    if min_prob is not None:
        hits = [hit for hit in hits if hit['Probab'] >= min_prob]
    if max_pvalue is not None:
        hits = [hit for hit in hits if hit['P-value'] <= max_pvalue]
    if max_evalue is not None:
        hits = [hit for hit in hits if hit['E-value'] <= max_evalue]
    if min_fragment_length is not None:
        hits = [hit for hit in hits if frag_size(hit) >= min_fragment_length]

    # read the original protein file, used to run HHsearch
    p = Protein.read(fullsequence_fp, seq_num=1)
    query_id = p.metadata['id']
    query_desc = p.metadata['description']

    results = {'match': [], 'non_match': []}
    # select non overlapping positive hits
    subseqs_pos = select_hits(hits, e_value_threshold=999999)

    for hit in subseqs_pos:
        _id = get_q_id(hit)
        match_id = hit['Hit'].split()[0]
        header = "%s %s %s" % (correct_header_positions(
            query_id,
            hit['alignment'][_id]['start'],
            hit['alignment'][_id]['end']), '# %s' % match_id, query_desc)
        seq = hit['alignment'][_id]['sequence'].replace('-', '')
        results['match'].append((header, seq, hit['alignment'][_id]['start']))

    # collect gaps between positive hits
    subseqs_neg = report_uncovered_subsequences(subseqs_pos, str(p),
                                                min_fragment_length)
    for hit in subseqs_neg:
        header = "%s %s" % (correct_header_positions(
            query_id,
            hit['start'],
            hit['end']), query_desc)
        seq = hit['sequence']
        results['non_match'].append((header, seq, hit['start']))

    # write sub-sequences to a multiple fasta file, sequences are un-wrapped
    try:
        # sort by start position
        for type_ in results:
            results[type_] = sorted(results[type_], key=lambda x: x[2])

        if subsequences_fp is not None:
            for type_ in results:
                f = open('%s.%s' % (subsequences_fp, type_), 'w')
                for res in results[type_]:
                    f.write(">%s\n%s\n" % res[:2])
                f.close()

        # removing the start position component from all subsequences
        return {type_: list(map(lambda x: x[:2], results[type_]))
                for type_ in results}
    except IOError:
        raise IOError('Cannot write to file "%s"' % subsequences_fp)


def pretty_output(mask_out):
    for key in sorted(mask_out.keys()):
        print(key)
        for i in range(len(mask_out[key])):
            print('\t>%s\n\t%s' % (mask_out[key][i][0], mask_out[key][i][1]))


# RUN FROM COMMAND LINE
@click.command()
@click.option('--subseq_fp', '-o', default=None,
              type=click.Path(),
              help='Root of output file\n'
              'Filepath to which sub-sequences are written as a multiple \
              fasta files. Each sequence makes up one header and one sequence \
              file, i.e. sequences are not wrapped.'
              'Two files will be produced, suffixed by ''.match'' and \
              ''.non_match''. The first holds sub-sequences of hits, the \
              second holds the none-hit covered subsequences.')
@click.option('--prob', '-p', default=None, type=float,
              help='Minimum HHsuite probability value')
@click.option('--e_val', '-e', default=None, type=float,
              help='Maximum E-value')
@click.option('--frag_len', '-l', default=0, type=float,
              help='Minimum fragment length')
@click.option('--p_val', '-v', default=None, type=float,
              help='Maximum P-value')
@click.argument('hh_fp', nargs=1, type=click.Path(exists=True))
@click.argument('fullseq_fp', nargs=1, type=click.Path(exists=True))
def _split_search(hh_fp, fullseq_fp, subseq_fp,
                  prob,
                  p_val,
                  e_val,
                  frag_len):

    _out = mask_sequence(hh_fp,
                         fullseq_fp,
                         subseq_fp,
                         prob,
                         p_val,
                         e_val,
                         frag_len)

    if subseq_fp is None:
        pretty_output(_out)


if __name__ == "__main__":
    _split_search()
