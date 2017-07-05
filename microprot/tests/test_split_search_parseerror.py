from unittest import TestCase, main

from skbio.util import get_data_path

from microprot.scripts.split_search import (mask_sequence,
                                            parse_pdb_match,
                                            _parse_hit_block)


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
        mask_sequence(self.file_hhsearch1,
                      self.file_fasta1,
                      min_prob=95.0,
                      min_fragment_length=40)

        parse_pdb_match(self.file_hhsearch1)

    def test_split_search_parseerror_2(self):
        mask_sequence(self.file_hhsearch2,
                      self.file_fasta2,
                      min_prob=95.0,
                      min_fragment_length=40)

        parse_pdb_match(self.file_hhsearch2)

    def test_split_search_parseerror_block(self):
        f = open(
            get_data_path('test_split_search/parsealignment_fail_example.txt'),
            'r')
        errorblock = "".join(f.readlines())
        f.close()
        _parse_hit_block(errorblock)


if __name__ == '__main__':
    main()
