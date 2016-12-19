from unittest import TestCase, main
from skbio.util import get_data_path
import re

from microprot.scripts.contacts import _topN_contacts
from microprot.scripts.contacts import (read_PDB_coordinates,
                                        read_contact_predictions,
                                        find_PDB_contacts,
                                        contact_precision,
                                        _contacts)


class ContactsTests(TestCase):
    def setUp(self):
        self.jbe_con = get_data_path('test_contacts/1jbeA.psicov')
        self.jbe_pdb = get_data_path('test_contacts/1jbeA_clean.pdb')

        self.qjp_con = get_data_path('test_contacts/1qjpA.psicov')
        self.qjp_pdb = get_data_path('test_contacts/1qjpA_clean.pdb')

        self.real_n_con_qjp = {'lr': 6441,
                               'sr': 2337,
                               'all': 8778}

        self.real_n_con_jbe = {'lr': 4742,
                               'sr': 2025,
                               'all': 6767}

        self.positive_params = [
            {'-t': 'all', '-l': 1},
            {'-t': 'all', '-l': 2},
            {'-t': 'all', '-l': 5},
            {'-t': 'all', '-l': 10},
            {'-t': 'lr', '-l': 1},
            {'-t': 'lr', '-l': 2},
            {'-t': 'lr', '-l': 10},
            {'-t': 'sr', '-l': 2},
            {'-t': 'sr', '-l': 10}]

    def test_read_PDB(self):
        for inp_fp, n_res in zip([self.jbe_pdb, self.qjp_pdb], [126, 137]):
            inp_fh = open(inp_fp, 'r')
            out = read_PDB_coordinates(inp_fh, 'CB')
            self.assertEqual(len(out), n_res)
            self.assertEqual(len(out[0]), 4)

    def test_read_contacts(self):
        for inp_fp, real_con in zip([self.jbe_con, self.qjp_con],
                                    [self.real_n_con_jbe, self.real_n_con_qjp]):
            for contype in ['all', 'sr', 'lr']:
                inp_fh = open(inp_fp, 'r')
                out_aa, out_ppv = read_contact_predictions(inp_fh,
                                                           contype=contype)
                self.assertEqual(len(out_aa), len(out_ppv))
                self.assertEqual(len(out_aa), real_con[contype])
            inp_fh.close()

    def test_find_contacts(self):
        for inp_fp in [self.jbe_pdb, self.qjp_pdb]:
            inp_fh = open(inp_fp, 'r')
            for contype in ['all', 'lr']:
                con_file = re.sub('_clean.pdb', '.contacts_', inp_fp)+contype
                for _c, _l in zip(find_PDB_contacts(read_PDB_coordinates(inp_fh),
                                                    con_type=contype),
                                  open(con_file).readlines()):
                    _i = _l.split()
                    self.assertEqual(_c[0], _i[0])
                    self.assertEqual(_c[1], _i[1])
                    self.assertEqual('%.4f' % _c[2], _i[2])
            inp_fh.close()

    def test_contact_precision(self):
        for con_fp, pdb_fp in zip([self.jbe_con, self.qjp_con],
                                  [self.jbe_pdb, self.qjp_pdb]):
            ref_fp = re.sub('.psicov', '.precision', con_fp)
            ref_lines = open(ref_fp, 'r').readlines()
            coords = read_PDB_coordinates(open(pdb_fp, 'r'))
            for params, ref in zip(self.positive_params, ref_lines):
                topN = _topN_contacts(coords, params['-l'])
                cons, ppvs = read_contact_predictions(open(con_fp, 'r'),
                                                      topX=topN,
                                                      contype=params['-t'])
                out = contact_precision(coords, cons, ppvs)

                self.assertEqual(ref.split()[4], '%.4f' % out)

if __name__ == '__main__':
    main()
