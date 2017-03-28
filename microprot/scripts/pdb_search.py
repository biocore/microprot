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
        content = lines[idx+1].split()[3]
        startCol = lines[idx+1].find(content)
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
            # check if there is an overlap to previous results
            noOverlap = True
            for good in good_hits:
                if is_overlapping((hit['alignment']['Q Consensus']['start'],
                                   hit['alignment']['Q Consensus']['end']),
                                  (good['alignment']['Q Consensus']['start'],
                                   good['alignment']['Q Consensus']['end'])):
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
        q_keys = [k for k in hit['alignment'].keys()
                  if k.startswith('Q') & (k != 'Q Consensus')]

        # compose a new dict as report for this hit.
        report.append({'pdb_id': hit['Hit'].split()[0],
                       'covered_sequence': hit['alignment'][q_keys[0]]
                       ['sequence'].replace('-', ''),
                       'start': hit['alignment'][q_keys[0]]['start'],
                       'end': hit['alignment'][q_keys[0]]['end']})
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
        q_keys = [k for k in hit['alignment'].keys()
                  if k.startswith('Q') & (k != 'Q Consensus')]
        covered.append((hit['alignment'][q_keys[0]]['start'],
                        hit['alignment'][q_keys[0]]['end']))

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


def mask_sequence(hhsuite_fp, fullsequence_fp, subsequences_fp=None,
                  min_prob=None, max_pvalue=None, max_evalue=None,
                  min_fragment_length=None):
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
        Default: None, i.e. no filtering on fragment length.

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
        def frag_size(hit):
            return (hit['alignment']['Q Consensus']['end'] -
                    hit['alignment']['Q Consensus']['start'] +
                    1)
        hits = [hit for hit in hits if frag_size(hit) >= min_fragment_length]

    # read the original protein file, used to run HHsearch
    p = Protein.read(fullsequence_fp, seq_num=1)
    queryname = p.metadata['id'] + " " + p.metadata['description']

    results = []
    # select non overlapping positive hits
    subseqs_pos = select_hits(hits, e_value_threshold=999999)
    for hit in subseqs_pos:
        header = "%s_%i-%i" % (queryname,
                               hit['alignment']['Q Consensus']['start'],
                               hit['alignment']['Q Consensus']['end'])
        seq = hit['alignment']['Q Consensus']['sequence'].replace('-', '')
        results.append((header, seq, hit['alignment']['Q Consensus']['start']))

    # collect gaps between positive hits
    subseqs_neg = report_uncovered_subsequences(subseqs_pos, str(p),
                                                min_subseq_len=0)
    for hit in subseqs_neg:
        header = "%s_%i-%i" % (queryname,
                               hit['start'],
                               hit['end'])
        seq = hit['sequence']
        results.append((header, seq, hit['start']))

    # write sub-sequences to a multiple fasta file, sequences are un-wrapped
    try:
        # sort by start position
        results = sorted(results, key=lambda x: x[2])

        if subsequences_fp is not None:
            f = open(subsequences_fp, 'w')
            for res in results:
                f.write(">%s\n%s\n" % res[:2])
            f.close()

        return list(map(lambda x: x[:2], results))
    except IOError:
        raise IOError('Cannot write to file "%s"' % subsequences_fp)
