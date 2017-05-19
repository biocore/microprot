from unittest import TestCase, main
from click.testing import CliRunner
from shutil import rmtree
from os import remove
from os.path import join, basename
from tempfile import mkdtemp
from skbio.util import get_data_path
from skbio import Sequence
from glob import glob

from microprot.scripts.process_fasta import (extract_sequences,
                                             write_sequences,
                                             read_representatives,
                                             split_fasta,
                                             _processing)


class ProcessingTests(TestCase):
    def setUp(self):
        # temporary working directory
        self.working_dir = mkdtemp()

        # test data files
        dir = 'test_process_fasta'
        self.input_faa = get_data_path(join(dir, 'input.faa'))
        self.represent = get_data_path(join(dir, 'represent.xml'))
        self.pdb_seqres = get_data_path(join(dir, 'pdb_seqres.txt'))

        self.split_1 = get_data_path(join(dir, '1K5N_B.fasta'))
        self.split_2 = get_data_path(join(dir, '2VB1_A.fasta'))
        self.split_3 = get_data_path(join(dir, '3J4F_A.fasta'))

        # test protein sequences
        self.seqs = []
        x = Sequence('PIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTML'
                     'NTVGGHQAAMQMLKETINEEAAEWDRLHPVHAGPIEPGQMREPRGSDIAGTTSTLQ'
                     'EQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRF'
                     'YKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGP'
                     'GHKARVL')
        x.metadata['id'] = '3J4F_A'
        x.metadata['description'] = ('Chain A, Structure Of Hiv-1 Capsid Prote'
                                     'in By Cryo-em')
        self.seqs.append(x)
        x = Sequence('MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLS'
                     'FSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM')
        x.metadata['id'] = '1K5N_B'
        x.metadata['description'] = ('Chain B, Hla-B2709 Bound To Nona-Peptide'
                                     ' M9')
        self.seqs.append(x)
        x = Sequence('KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGIL'
                     'QINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWR'
                     'NRCKGTDVQAWIRGCRL')
        x.metadata['id'] = '2VB1_A'
        x.metadata['description'] = ('Chain A, Hewl At 0.65 Angstrom Resolutio'
                                     'n')
        self.seqs.append(x)

    def test_extract_sequences(self):
        # extract all sequences
        obs = extract_sequences(self.input_faa)
        exp = self.seqs
        self.assertListEqual(obs, exp)

        # specify protein index
        obs = extract_sequences(self.input_faa, identifiers=2)
        exp = [self.seqs[1]]
        self.assertListEqual(obs, exp)

        # specify protein ID
        obs = extract_sequences(self.input_faa, identifiers='1K5N_B')
        exp = [self.seqs[1]]
        self.assertListEqual(obs, exp)

        # specify protein indexes using a Python list of int
        obs = extract_sequences(self.input_faa, identifiers=[1, 3])
        exp = [self.seqs[0], self.seqs[2]]
        self.assertListEqual(obs, exp)

        # specify protein indexes using a Python list of str
        obs = extract_sequences(self.input_faa, identifiers=['2', '3'])
        exp = [self.seqs[1], self.seqs[2]]
        self.assertListEqual(obs, exp)

        # specify protein indexes using a comma-separated list
        obs = extract_sequences(self.input_faa, identifiers='1,3')
        exp = [self.seqs[0], self.seqs[2]]
        self.assertListEqual(obs, exp)

        # specify protein entries using an external file
        listfile = join(self.working_dir, 'list.txt')
        with open(listfile, 'w') as f:
            f.write('%s\n%s\n' % ('1K5N_B', '2VB1_A'))
        obs = extract_sequences(self.input_faa, identifiers=listfile)
        remove(listfile)
        exp = [self.seqs[1], self.seqs[2]]
        self.assertListEqual(obs, exp)

        # specify protein index range using a Python tuple (start, end)
        obs = extract_sequences(self.input_faa, identifiers=(2, 3))
        exp = [self.seqs[1], self.seqs[2]]
        self.assertListEqual(obs, exp)

        # raise error when tuple is not properly formatted
        err = 'Error: Index range must be a tuple of (start, end).'
        with self.assertRaises(ValueError, msg=err):
            extract_sequences(self.input_faa, identifiers=('hi', 'there'))
        with self.assertRaises(ValueError, msg=err):
            extract_sequences(self.input_faa, identifiers=(1, 2, 3))
        with self.assertRaises(ValueError, msg=err):
            extract_sequences(self.input_faa, identifiers=(3, 1))

        # specify protein index range using a str of "start..end"
        obs = extract_sequences(self.input_faa, identifiers='1..3')
        exp = [self.seqs[0], self.seqs[1], self.seqs[2]]
        self.assertListEqual(obs, exp)

        # raise error when "start..end" is not properly formatted
        err = 'Error: Index range must be formatted as "start..end".'
        with self.assertRaises(ValueError, msg=err):
            extract_sequences(self.input_faa, identifiers='not.an..id')
        with self.assertRaises(ValueError, msg=err):
            extract_sequences(self.input_faa, identifiers='1..2..3')
        with self.assertRaises(ValueError, msg=err):
            extract_sequences(self.input_faa, identifiers='3..1')

        # raise error when identifiers is of incorrect data type
        err = 'Error: Incorrect data type of identifiers.'
        with self.assertRaises(ValueError, msg=err):
            extract_sequences(self.input_faa, identifiers=1.23)
        with self.assertRaises(ValueError, msg=err):
            extract_sequences(self.input_faa, identifiers={'1K5N_B': 1})

    def test_write_sequences(self):
        seqs = extract_sequences(self.input_faa, identifiers='1')
        outfile = join(self.working_dir, 'output.faa')
        write_sequences(seqs, outfile)
        with open(outfile, 'r') as f:
            obs = f.read().splitlines()
        remove(outfile)
        exp = ['>%s %s' % (self.seqs[0].metadata['id'],
                           self.seqs[0].metadata['description']),
               str(self.seqs[0])]
        self.assertListEqual(obs, exp)

    def test_read_representatives(self):
        # read from local file
        obs = read_representatives(self.represent)
        exp = ['1k5n_B', '2vb1_A', '3j4f_A']
        self.assertListEqual(obs, exp)

        # read from server (test if not empty)
        obs = read_representatives('90')
        self.assertTrue(obs)

        # bad parameter
        err = ('Error: You must specify a local file path or the clustering '
               'level.')
        with self.assertRaisesRegex(ValueError, err):
            read_representatives('invalid_string')

    def test_split_fasta(self):
        # without prefix
        for prefix in [None, 'dupa']:
            for inp_faa in [extract_sequences(self.input_faa),
                            self.input_faa]:
                split_fasta(inp_faa, prefix=prefix,
                            outdir=self.working_dir)
                if prefix is not None:
                    exp_fnames = sorted(['%s_%s' % (prefix,
                                                    basename(self.split_1)),
                                         '%s_%s' % (prefix,
                                                    basename(self.split_2)),
                                         '%s_%s' % (prefix,
                                                    basename(self.split_3))])
                else:
                    exp_fnames = sorted([basename(self.split_1),
                                         basename(self.split_2),
                                         basename(self.split_3)])

                obs_paths = sorted(glob(self.working_dir+'/*.fasta'))
                obs_fnames = list(map(basename,
                                      obs_paths))
                self.assertListEqual(obs_fnames, exp_fnames)

                for obs_fp, exp_fp in zip([self.split_1, self.split_2,
                                          self.split_3], obs_paths):
                    obs = open(obs_fp, 'r').read()
                    exp = open(exp_fp, 'r').read()
                    self.assertEqual(obs, exp)

                # remove fastas
                for fasta in obs_paths:
                    remove(fasta)

    def test__processing(self):
        # get representative proteins
        params = ['--infile', self.input_faa,
                  '--outfile', join(self.working_dir, 'output.faa'),
                  '--identifiers', None,
                  '--represent', self.represent]
        res = CliRunner().invoke(_processing, params)
        self.assertEqual(res.exit_code, 0)
        exp = ('Number of representative proteins: 3\n'
               'Number of extracted proteins: 0\n'
               'Task completed.\n')
        self.assertEqual(res.output, exp)

        # extract proteins
        params = ['--infile', self.input_faa,
                  '--outfile', join(self.working_dir, 'output.faa'),
                  '--identifiers', '1,3',
                  '--represent', None]
        res = CliRunner().invoke(_processing, params)
        self.assertEqual(res.exit_code, 0)
        exp = ('Number of extracted proteins: 2\n'
               'Task completed.\n')
        self.assertEqual(res.output, exp)

        # TODO split fasta

    def tearDown(self):
        rmtree(self.working_dir)


if __name__ == '__main__':
    main()
