from unittest import TestCase, main

from skbio.util import get_data_path

from microprot.scripts.split_search import (mask_sequence)


class ParsersTests(TestCase):
    def setUp(self):
        self.fp_out = get_data_path('test_split_search/NC_000913.3_2.out')
        self.fp_seqs = get_data_path('test_split_search/NC_000913.3_2.fasta')

    def test_mask_sequence(self):
        match_header_2 = (
            'gi|556503834|ref|NC_000913.3|_2_464-815 # 1ebf_A # 337 # 2799 # '
            '1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spa'
            'cer=5-10bp;gc_cont=0.531')

        # test default behaviour
        obs = mask_sequence(self.fp_out, self.fp_seqs)
        self.assertEqual([m[0] for m in obs['match']], [
            ('gi|556503834|ref|NC_000913.3|_2_1-461 # 2j0w_A # 337 # 2799 # 1 '
             '# ID=1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_space'
             'r=5-10bp;gc_cont=0.531'), match_header_2])

        # restrict hits to satisfy at least 33% sequence identifty
        # --> first match will differ from default above
        obs = mask_sequence(self.fp_out, self.fp_seqs, min_identity=0.33)
        self.assertEqual([m[0] for m in obs['match']], [
            ('gi|556503834|ref|NC_000913.3|_2_1-462 # 3c1m_A # 337 # 2799 # 1 '
             '# ID=1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_space'
             'r=5-10bp;gc_cont=0.531'), match_header_2])


if __name__ == '__main__':
    main()
