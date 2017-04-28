# edge cases:
## input "runs out" before we reach the end of the pipeline
## input has only 1 domain
## sequence contains unknown AAs


from unittest import TestCase, main

from skbio.util import get_data_path

from microprot.scripts.pdb_search import is_overlapping, \
                                         _parse_hit_summary_line, \
                                         _parse_hit_block, parse_pdb_match, \
                                         select_hits, report_hits, \
                                         report_uncovered_subsequences


class ParsersTests(TestCase):
    def setUp(self):
        self.line_a = ('  3 2cdq_A Aspartokinase; aspartat 100.0 4.1E-58 '
                       '1.1E-62  453.8   0.0  458    1-473    28-496 (510)\n')
        self.true_a = {'No': 3,
                       'Hit': '2cdq_A Aspartokinase; aspartat',
                       'Prob': 100.0,
                       'E-value': 4.1e-58,
                       'P-value': 1.1e-62,
                       'Score': 453.8,
                       }
        self.line_b = ('    36 1gtm_A Glutamate dehydrogenase  98.8 2.9E-12 '
                       '7.9E-17  126.1   0.0  201  451-679   193-407 (419)\n')
        self.file_a = get_data_path('test_pdb_search/NC_000913.3_2.out')
        self.file_b = get_data_path('test_pdb_search/T0810-D1.fasta.out')
        self.true_subseqs = [{'sequence': 'TD', 'start': 462, 'end': 463},
                             {'sequence': 'WKLGV', 'start': 816, 'end': 820}]

    def test_is_overlapping(self):
        self.assertTrue(is_overlapping((20, 100), (30, 150)))
        self.assertTrue(is_overlapping((20, 100), (2, 70)))

        with self.assertRaises(ValueError) as ve:
            is_overlapping((100, 20), (1, 8))
        self.assertEqual(
            ('Left component of intA cannot be larger than its right component'
             ', i.e. the interval must be forward oriented.'),
            str(ve.exception)
        )
        with self.assertRaises(ValueError) as ve:
            is_overlapping((10, 20), (-9, -20))
        print(str(ve.exception))
        self.assertEqual(
            ('Left component of intB cannot be larger than its right component'
             ', i.e. the interval must be forward oriented.'),
            str(ve.exception)
        )

    def test__parse_hit_summary_line(self):
        self.assertEqual(self.true_a, _parse_hit_summary_line(self.line_a))
        self.assertEqual(self.true_b, _parse_hit_summary_line(self.line_b))
        with self.assertRaises(ValueError):
            _parse_hit_summary_line('no valid line')

    def test__parse_hit_block(self):
        self.assertEqual(self.true_block_a, _parse_hit_block(self.block_a))

    def test_parse_pdb_match(self):
        res = parse_pdb_match(self.file_b)
        self.assertEqual(len(res), 10)

        r0 = res[0]
        del r0['Cols']
        del r0['P-value']
        del r0['SS']
        self.assertEqual(self.true_block_a, r0)

        r9 = res[9]
        del r9['Cols']
        del r9['P-value']
        del r9['SS']
        self.assertEqual(self.true_block_b, r9)

        with self.assertRaises(IOError):
            parse_pdb_match('/does/not/exist')

    def test_select_hits(self):
        hits = parse_pdb_match(self.file_b)
        self.assertEqual(len(select_hits(hits, e_value_threshold=100)), 2)
        self.assertEqual(len(select_hits(hits)), 0)

        hits = parse_pdb_match(self.file_a)
        goodhits = select_hits(hits)
        self.assertEqual(len(goodhits), 2)
        self.assertEqual(goodhits[0]['No'], 1)
        self.assertEqual(goodhits[1]['No'], 5)

    def test_report_hits(self):
        hits = select_hits(parse_pdb_match(self.file_a))
        rep = report_hits(hits)
        self.assertEqual(self.true_a_hits_report, rep)

    def test_report_uncovered_subsequences(self):
        hits = select_hits(parse_pdb_match(self.file_a))
        subseqs = report_uncovered_subsequences(hits, self.query, 0)
        self.assertEqual(subseqs, self.true_subseqs)

if __name__ == '__main__':
    main()
