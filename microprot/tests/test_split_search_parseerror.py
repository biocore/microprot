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
        self.errorblock = """No 92
>1UYN_X NALP; AUTOTRANSPORTER, TRANSLOCATOR DOMAIN, MEMBRANE PROTEIN; HET: CXE, SO4; 2.6A {NEISSERIA MENINGITIDIS} SCOP: f.4.5.1
Probab=22.45  E-value=52  Score=25.20  Aligned_cols=66  Identities=12%  Similarity=0.115  Sum_probs=35.1  Template_Neff=11.100

Q GRAMNEG_T1D_31  114 LFEIEAVAGMGWLHYYVNG-------------DGDQNSWSTRFGLNFNFNLGETKAWTLGIKPAIVYDMQGNFPETKSRF  180 (275)
Q Consensus       114 lfeieavagmgwlhyyvng-------------dgdqnswstrfglnfnfnlgetkawtlgikpaivydmqgnfpetksrf  180 (275)
                      -+.|+..+++.+.+....+             ..+.+++..+.|+.+...+..  .|.+-.+-...++..+.-+....+|
T Consensus       170 ~~~i~P~~~l~~~~~~~~~~~e~g~~~~l~~~~~~~~~~~~~~G~~~~~~~~~--~~~~~~~~~~~~~~~~~~~~~~~~~  247 (308)
T 1UYN_X          170 TGDLTVEGGLRYDLLKQDAFAEKGSALGWSGNSLTEGTLVGLAGLKLSQPLSD--KAVLFATAGVERDLNGRDYTVTGGF  247 (308)
T ss_dssp             CCEEEEEEEEEEEEEEECCEEEESCTTCCEECCEEEEEEEEEEEEEEEEECSS--SEEEEEEEEEEEESSCCCCCCCC--
T ss_pred             ceEEEEEEEEEEEEeeccceeeeccccceEeeecccceEEEEeEEEEEeeccc--eEEEEEEEEEEEecCCCceEEeeEe
Confidence            3556666666665543322             123456667777777766654  4555555556666655544444444


Q GRAMNEG_T1D_31  181 N  181 (275)
Q Consensus       181 n  181 (275)
                      .
T Consensus       248 ~  248 (308)
T 1UYN_X          248 T  248 (308)
T ss_dssp             -
T ss_pred             e
Confidence            4


"""

    def test_split_search_parseerror_1(self):
        with self.assertRaisesRegex(ValueError,
                                    "with base 10: 'EG_T1D_51'"):
            mask_sequence(self.file_hhsearch1,
                          self.file_fasta1,
                          min_prob=95.0,
                          min_fragment_length=40)

        with self.assertRaisesRegex(ValueError,
                                    "with base 10: 'EG_T1D_51'"):
            parse_pdb_match(self.file_hhsearch1)

    def test_split_search_parseerror_2(self):
        with self.assertRaisesRegex(ValueError,
                                    "with base 10: 'MNEG_T1D_31'"):
            mask_sequence(self.file_hhsearch2,
                          self.file_fasta2,
                          min_prob=95.0,
                          min_fragment_length=40)

        with self.assertRaisesRegex(ValueError,
                                    "with base 10: 'MNEG_T1D_31'"):
            parse_pdb_match(self.file_hhsearch2)

    def test_split_search_parseerror_block(self):
        with self.assertRaisesRegex(ValueError,
                                    "with base 10: 'MNEG_T1D_31'"):
            _parse_hit_block(self.errorblock)


if __name__ == '__main__':
    main()
