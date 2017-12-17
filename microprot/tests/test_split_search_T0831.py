from unittest import TestCase, main

from skbio.util import get_data_path

from microprot.scripts.split_search import (is_overlapping,
                                            _parse_hit_summary_line,
                                            _parse_hit_block, parse_pdb_match,
                                            select_hits, report_hits,
                                            report_uncovered_subsequences,
                                            frag_size, get_q_id, mask_sequence)


class ParsersTests(TestCase):
    def setUp(self):
        self.fp_out = get_data_path('test_split_search/T0831.out')
        self.fp_seqs = get_data_path('test_split_search/T0831.fna')

    def test_mask_sequence(self):
        obs = mask_sequence(self.fp_out, self.fp_seqs, subsequences_fp='kurt_',
                            min_prob=95.0, max_evalue=0.1,
                            min_fragment_length=40)
        print(obs)

if __name__ == '__main__':
    main()
