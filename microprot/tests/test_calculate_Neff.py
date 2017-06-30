from unittest import TestCase, main
from click.testing import CliRunner
from shutil import rmtree
from os.path import join
from tempfile import mkdtemp
from skbio.util import get_data_path
from skbio import DistanceMatrix

from microprot.scripts.calculate_Neff import (parse_msa_file,
                                              hamming_distance_matrix,
                                              cluster_sequences,
                                              effective_family_size,
                                              _calculate_Neff)


class ProcessingTests(TestCase):
    def setUp(self):
        # temporary working directory
        self.working_dir = mkdtemp()
        self.output_fp = join(self.working_dir, 'output.txt')

        # test data files
        dir = 'test_calculate_Neff'
        self.input_a3m_fp = get_data_path(join(dir, '2phyA.a3m'))
        self.input_single_a3m_fp = get_data_path(join(dir, 'single.a3m'))
        self.plain_msa_fp = get_data_path(join(dir, '2phyA.aln'))
        self.hamming_dm_fp = get_data_path(join(dir, '2phyA.hdm'))
        self.clusters_fp = get_data_path(join(dir, '2phyA.c80'))

    def test_parse_msa_file(self):
        obs = [str(x) for x in parse_msa_file(self.input_a3m_fp)]
        with open(self.plain_msa_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)

    def test_hamming_distance_matrix(self):
        msa = parse_msa_file(self.input_a3m_fp)
        obs = hamming_distance_matrix(msa)
        exp = DistanceMatrix.read(self.hamming_dm_fp)
        self.assertEqual(obs, exp)

    def test_cluster_sequences(self):
        msa = parse_msa_file(self.input_a3m_fp)
        hdm = hamming_distance_matrix(msa)
        obs = cluster_sequences(hdm, 80)
        with open(self.clusters_fp, 'r') as f:
            exp = [int(x.split('\t')[1]) for x in f.read().splitlines()]
        self.assertListEqual(obs, exp)
        hdm = DistanceMatrix([[0]])
        msg = ('The number of observations cannot be determined on an empty '
               'distance matrix.')
        with self.assertRaisesRegex(ValueError, msg):
            cluster_sequences(hdm, 80)

    def test_effective_family_size(self):
        msa = parse_msa_file(self.input_a3m_fp)
        hdm = hamming_distance_matrix(msa)
        clu = cluster_sequences(hdm, 80)
        obs = effective_family_size(clu, msa.shape[1])
        exp = 6.6187612134
        self.assertAlmostEqual(obs, exp)
        obs = effective_family_size(clu, msa.shape[1], Nclu=True)
        exp = 74
        self.assertEqual(obs, exp)

    def test__calculate_Neff(self):
        params = ['--infile', self.input_a3m_fp,
                  '--outfile', self.output_fp,
                  '--cutoff', 80]
        res = CliRunner().invoke(_calculate_Neff, params)
        self.assertIn('Effective family size at 80% identity: 6.619.',
                      res.output)
        with open(self.output_fp, 'r') as f:
            obs = float(f.readline().strip())
        exp = 6.6187612134
        self.assertAlmostEqual(obs, exp)
        self.assertEqual(res.exit_code, 0)
        self.assertIn('Task completed.', res.output)
        params = ['--infile', self.input_single_a3m_fp]
        res = CliRunner().invoke(_calculate_Neff, params)
        self.assertIn('Effective family size at 100% identity: 0.000.',
                      res.output)

    def tearDown(self):
        rmtree(self.working_dir)


if __name__ == '__main__':
    main()
