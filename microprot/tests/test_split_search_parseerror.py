from unittest import TestCase, main

from skbio.util import get_data_path

from microprot.scripts.split_search import (mask_sequence,
                                            correct_header_positions)


class ParsersTests(TestCase):
    def setUp(self):
        self.file_hhsearch1 = get_data_path(
            'test_split_search/GRAMNEG_T1D_5168.out')
        self.file_fasta1 = get_data_path(
            'test_split_search/GRAMNEG_T1D_5168.fasta')
        self.file_hhsearch2 = get_data_path(
            'test_split_search/GRAMNEG_T1D_3144_1-275.out')
        self.file_fasta2 = get_data_path(
            'test_split_search/GRAMNEG_T1D_3144_1-275.fasta')

    def test_split_search_parseerror_1(self):
        with self.assertRaisesRegex(ValueError,
                                    "with base 10: 'EG_T1D_51'"):
            mask_sequence(self.file_hhsearch1,
                          self.file_fasta1,
                          min_prob=95.0,
                          min_fragment_length=40)

    def test_split_search_parseerror_2(self):
        with self.assertRaisesRegex(ValueError,
                                    "with base 10: 'MNEG_T1D_31'"):
            mask_sequence(self.file_hhsearch2,
                          self.file_fasta2,
                          min_prob=95.0,
                          min_fragment_length=40)


if __name__ == '__main__':
    main()
