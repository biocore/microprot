from unittest import TestCase, main
import tempfile
import os
import shutil

from skbio.util import get_data_path

from microprot.scripts.split_search import (mask_sequence, parse_pdb_match)


class ParsersTests(TestCase):
    def setUp(self):
        self.file_missed = get_data_path(
            'test_split_search/NZ_GG666849.1_2_251-330.out')
        self.file_query = get_data_path(
            'test_split_search/NZ_GG666849.1_2_251-330.fasta')
        self.exp_hits = {
            'non_match':
            [(('NZ_GG666849.1_2_251-330_1-80 # 798 # 2885 # -1 # ID=1_2;partia'
               'l=00;start_type=TTG;rbs_motif=AGxAGG/AGGxGG;rbs_spacer=5-10bp;'
               'gc_cont=0.499'),
              ('IGIQGDTYSEDEDYPELPRTANGRLSSYILVNHKEQVHVYNQIATKLGLQKESGEVVMLPSQ'
               'FINRFSLRNEHGRGIPDQ'))],
            'match': []}

    def test_parsing(self):
        obs_res = parse_pdb_match(self.file_missed)
        self.assertEqual(len(obs_res), 12)

    def test_mask_sequence(self):
        obs_res = mask_sequence(self.file_missed,
                                self.file_query,
                                None,
                                max_evalue=0.1,
                                min_fragment_length=40)

        # check that hit selection works correct
        self.assertEqual(obs_res, self.exp_hits)

    def test_outdir(self):
        outdir = tempfile.mkdtemp(prefix='splitseq_')

        mask_sequence(self.file_missed,
                      self.file_query,
                      outdir,
                      max_evalue=0.1,
                      min_fragment_length=40)

        exp_filecontents = {
            'non_match':
            [('>NZ_GG666849.1_2_251-330_1-80 # 798 # 2885 # -1 # ID=1_2;parti'
              'al=00;start_type=TTG;rbs_motif=AGxAGG/AGGxGG;rbs_spacer=5-10bp'
              ';gc_cont=0.499\n'),
             ('IGIQGDTYSEDEDYPELPRTANGRLSSYILVNHKEQVHVYNQIATKLGLQKESGEVVMLPSQ'
              'FINRFSLRNEHGRGIPDQ\n')],
            'match': []}
        obs_filecontents = dict()
        for type_ in ('match', 'non_match'):
            filename = outdir+'.'+type_
            f = open(filename, 'r')
            obs_filecontents[type_] = f.readlines()
            f.close()
            os.remove(filename)
        shutil.rmtree(outdir)

        self.assertEqual(obs_filecontents, exp_filecontents)


if __name__ == '__main__':
    main()
