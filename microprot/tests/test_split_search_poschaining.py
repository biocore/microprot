from unittest import TestCase, main

from skbio.util import get_data_path

from microprot.scripts.split_search import (mask_sequence,
                                            correct_header_positions)


class ParsersTests(TestCase):
    def setUp(self):
        self.file_out1 = get_data_path(
            'test_split_search/poschaining_1.out')
        self.file_fasta1 = get_data_path(
            'test_split_search/poschaining_1.fasta')
        self.pos1 = [('match', '223', '250'),
                     ('match', '331', '380'),
                     ('non_match', '1', '222'),
                     ('non_match', '251', '330'),
                     ('non_match', '381', '696')]

        self.file_out2 = get_data_path(
            'test_split_search/poschaining_2.out')
        self.file_fasta2 = get_data_path(
            'test_split_search/poschaining_2.fasta')
        self.pos2 = [('non_match', '251', '330')]

    def test_level1(self):
        obs = mask_sequence(self.file_out1, self.file_fasta1, max_pvalue=0.95,
                            min_fragment_length=10)

        positions = []
        for type_ in sorted(obs.keys()):
            for header in sorted(obs[type_]):
                start, stop = \
                    header[0].split(' # ')[0].split('_')[-1].split('-')
                positions.append((type_, start, stop))

        self.assertEqual(positions, self.pos1)

    def test_level2(self):
        obs = mask_sequence(self.file_out2, self.file_fasta2, max_evalue=0.1,
                            min_fragment_length=40)

        positions = []
        for type_ in sorted(obs.keys()):
            for header in sorted(obs[type_]):
                start, stop = \
                    header[0].split(' # ')[0].split('_')[-1].split('-')
                positions.append((type_, start, stop))

        self.assertEqual(positions, self.pos2)

    def test_updating(self):
        prefix = 'test_1_'
        self.assertEqual(
            correct_header_positions('%s251-330' % prefix, 30, 74),
            '%s280-324' % prefix)

        with self.assertRaises(ValueError) as ve:
            correct_header_positions('test_1_251-330', 300, 374)
        self.assertEqual(str(ve.exception),
                         "New interval out of old interval borders!")

        with self.assertRaises(ValueError) as ve:
            correct_header_positions('test_1_251-330', 30, 174)
        self.assertEqual(str(ve.exception),
                         "Intervals are not nested!")

        self.assertEqual(correct_header_positions('test_1_251330', 30, 74),
                         "test_1_251330_30-74")


if __name__ == '__main__':
    main()
