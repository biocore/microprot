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

        self.assertEqual(
            obs['match'][0][1],
            ('TMEELLTSLQKKCGTECEEAHRQLVCALNGLAGIHIIKGEYALAAELYREVLRSSEEHKGKLK'
             'TDSLQRLHATHNLMELLIARHPGIPPTLRDGRLEEEAKQLREHYMSKCNTEVAEAQQALYPVQ'
             'QTIHELQRKIHSNSPWWLNVIHRAIEFTIDEELVQRVRNEITSNYKQQTGKLSMSEKFRDCRG'
             'LQFLLTTQMEELNKCQKLVREAVKNLEGPPSRNVIESATVCHLRPARLPLNCCVFCKADELFT'
             'EYESKLFSNTVKGQTAIFEEMIEDEEGLVDDRAPTTTRGLWAISETERSMKAILSFAKSHRFD'
             'VEFVDEGSTSMDLFEAWKKEYKLLHEYWMALRNRVSAVDELAMATERLRVRDPREPKPNPPVL'
             'HIIEPHEVEQNRIKLLNDKAVATSQLQKKLGQLLYLTNLEK'))

        exp_0 = {
            'Probab': 100.0,
            'Template_Neff': 8.5,
            'P-value': 2.8e-85,
            'Similarity': 1.445,
            'Sum_probs': 363.5,
            'Score': 555.49,
            'Cols': 419,
            'No': 1,
            'Identities': 1.0,
            'SS': 0.0,
            'alignment': {
                'Q T0831': {
                    'start': 1,
                    'end': 419,
                    'sequence':
                    ('TMEELLTSLQKKCGTECEEAHRQLVCALNGLAGIHIIKGEYALAAELYREVLRSSE'
                     'EHKGKLKTDSLQRLHATHNLMELLIARHPGIPPTLRDGRLEEEAKQLREHYMSKCN'
                     'TEVAEAQQALYPVQQTIHELQRKIHSNSPWWLNVIHRAIEFTIDEELVQRVRNEIT'
                     'SNYKQQTGKLSMSEKFRDCRGLQFLLTTQMEELNKCQKLVREAVKNLEGPPSRNVI'
                     'ESATVCHLRPARLPLNCCVFCKADELFTEYESKLFSNTVKGQTAIFEEMIEDEEGL'
                     'VDDRAPTTTRGLWAISETERSMKAILSFAKSHRFDVEFVDEGSTSMDLFEAWKKEY'
                     'KLLHEYWMALRNRVSAVDELAMATERLRVRDPREPKPNPPVLHIIEPHEVEQNRIK'
                     'LLNDKAVATSQLQKKLGQLLYLTNLEK'),
                    'totallen': 419},
                'Q Consensus': {
                    'start': 1,
                    'end': 419,
                    'sequence':
                    ('tmeelltslqkkcgteceeahrqlvcalnglagihiikgeyalaaelyrevlrsse'
                     'ehkgklktdslqrlhathnlmelliarhpgipptlrdgrleeeakqlrehymskcn'
                     'tevaeaqqalypvqqtihelqrkihsnspwwlnvihraieftideelvqrvrneit'
                     'snykqqtgklsmsekfrdcrglqfllttqmeelnkcqklvreavknlegppsrnvi'
                     'esatvchlrparlplnccvfckadelfteyesklfsntvkgqtaifeemiedeegl'
                     'vddraptttrglwaisetersmkailsfakshrfdvefvdegstsmdlfeawkkey'
                     'kllheywmalrnrvsavdelamaterlrvrdprepkpnppvlhiiepheveqnrik'
                     'llndkavatsqlqkklgqllyltnlek'),
                    'totallen': 419},
                'column score': {
                    'sequence':
                    ('||+|++..|-++|-+|||+++|++|.++|||||||||+|+|..|+++||+||+..+'
                     '++++++++|+||++|+.|||.+++...+||+||+++|..+.+++.+++..|++++.'
                     '..+..|++.+.++.+.+++++.+.++.++||+.+++.+++..++..++++|+++++'
                     '.+|.+..|..++..+|++.+||.+.+++.+++|.++.+-+.+++++|++||..+++'
                     '+++..||++|.+-+...|.+|++++.|..||+.||+.+.+|.+..+++++++++|.'
                     '.++.....++|.|+.|+.|+.+|.|++|++++.|+.+++.+|..-++++++|||||'
                     '+.++.+|++.+..++|.|||.|++-|+|.++|.++.|+||..++|.|+++++.+.+'
                     '+.+++.++...|++++|||.||.||.|')},
                'T Consensus': {
                    'start': 2,
                    'end': 420,
                    'sequence':
                    ('tmeell~~Li~k~~~eceea~R~~v~~~NgLAgl~~l~~~~~~A~~~YrevL~~~~'
                     '~~~~~~~~D~Lq~iH~l~NL~~~l~~~~~~~~~~~~~~~l~~~~~~l~~~Yl~~~~'
                     '~~~~~a~~~~~~~~~~~~~~~~~~~~~~~Ww~~~l~~~~~~~~~~~l~~~i~~~l~'
                     '~~~~~~~~~~~~~~~~~s~~gL~~~l~~~l~~L~~~R~~l~~~l~~L~~~~~~~~v'
                     '~~~~~Ch~~~~~~~~~~C~~C~~~~~~~~yE~~Lf~~~~~~~~~~~~~~~~~~~~~'
                     '~~~~~~~~~~g~~~~S~~e~~lk~i~~~~r~~~~~~~~~~~~~~hl~~le~~rkEf'
                     '~~~r~lw~~~~~~l~a~DEL~ma~~Rlrl~~~~e~~~~~~~~~~i~~~ev~~~~~~'
                     '~~~e~~~a~~~l~r~~gqLrYL~nL~k'),
                    'totallen': 420},
                'T 4QN1_A': {
                    'start': 2,
                    'end': 420,
                    'sequence':
                    ('TMEELLTSLQKKCGTECEEAHRQLVCALNGLAGIHIIKGEYALAAELYREVLRSSE'
                     'EHKGKLKTDSLQRLHATHNLMELLIARHPGIPPTLRDGRLEEEAKQLREHYMSKCN'
                     'TEVAEAQQALYPVQQTIHELQRKIHSNSPWWLNVIHRAIEFTIDEELVQRVRNEIT'
                     'SNYKQQTGKLSMSEKFRDCRGLQFLLTTQMEELNKCQKLVREAVKNLEGPPSRNVI'
                     'ESATVCHLRPARLPLNCCVFCKADELFTEYESKLFSNTVKGQTAIFEEMIEDEEGL'
                     'VDDRAPTTTRGLWAISETERSMKAILSFAKSHRFDVEFVDEGSTSMDLFEAWKKEY'
                     'KLLHEYWMALRNRVSAVDELAMATERLRVRDPREPKPNPPVLHIIEPHEVEQNRIK'
                     'LLNDKAVATSQLQKKLGQLLYLTNLEK'),
                    'totallen': 420},
                'T ss_dssp': {
                    'sequence':
                    ('HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHTCHHHHHHHHHHHHHHHH'
                     'HTTTTCCCCHHHHHHHHHHHHHCCCCCTTSSCCCTTTTTHHHHHHHHHHHHHHHHH'
                     'HHHHHHHHTTHHHHHHHHHHHHSSCSSSCHHHHHHHHHHHTTCHHHHHHHHHHHHC'
                     'CC----------GGGCSSHHHHHHHHHHHHHHHHHHHHHHHHHHHTTCSSCCHHHH'
                     'HHHCCCCCSCSSSCCCCSHHHHHHHHHHHHHHHHBCCC------------------'
                     '-----------CCSBCHHHHHHHHHHHHHHHTTCCHHHHHHHHHHHHHHHHHHHHH'
                     'HHHHHHHHHHHHHHHHHHHHHHHHCCCEECCC---------CCEECTTCHHHHHHH'
                     'HHHHHHHHHHHHHHHHHHHHHHHTTCC')},
                'T ss_pred': {
                    'sequence':
                    ('CHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhCCHHHHHHHHHHHHHHHH'
                     'HhhcCCccchHHHHHHHhhHHHHHHhcCCCCCCCcchhHHHHHHHHHHHHHHHHHH'
                     'HHHHHHHHHHHHHHHHHHHHHHhhccCCcHHHHHHHHHHHCCCcHHHHHHHHHHHH'
                     'hhcccccCCcccccccccHHHHHHHHHHHHHHHHHHHHHHHHHHHhhcCCCcHHHH'
                     'HHhhcCCCCCCCCCCCCCCccccHHHHHHHHHHHhhcccCCCccchHhhhhccccc'
                     'cccCCCcccCCcccccHHHHHHHHHHHHHHhcCCCHHHHHHHHHHHHHHHHHHHHH'
                     'HHHHHHHHHHHHHHHHHHHHHHchhhheeCCCCCCCCCCCcccccCHHHHHHHHHH'
                     'HHHHHHHHHHHHHHHHHHHHHHhcccC')},
                'Confidence': {
                    'sequence':
                    ('79999999999999999999999999999999999999999999999999999999'
                     '99999999999999999999999999999999999999999999999999999999'
                     '99999999999999999999988888999999999998877899999999999999'
                     '99987777778999999999999999999999999999999999999999999999'
                     '99999999999877789999999999999999999999999999999999999999'
                     '99888889999999999999999999999999999999999999999999999999'
                     '99999999999999999999999999999999999999999999999999999999'
                     '999999999999999999999999976')},
            },
            'Aligned_cols': 419,
            'E-value': 1.6e-80,
            'Hit': ('4QN1_A E3 ubiquitin-protein ligase SHPRH; SHPRH, E3 ligas'
                    'e, RING, Ubiquitin; 2.48A {Homo sapiens}')
            }
        obs = parse_pdb_match(self.fp_out)
        for k in obs[0].keys():
            if type(obs[0][k]) == dict():
                self.assertCountEqual(obs[0][k], exp_0[k])
            else:
                self.assertEqual(obs[0][k], exp_0[k])

        with open(get_data_path('test_split_search/T0831_block0.out'),
                  'r') as f:
            block = "".join(f.readlines())
        obs = _parse_hit_block(block)
        for k in obs.keys():
            if type(obs[k]) == dict():
                self.assertCountEqual(obs[k], exp_0[k])
            else:
                self.assertEqual(obs[k], exp_0[k])


if __name__ == '__main__':
    main()
