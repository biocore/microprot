from unittest import TestCase, main

from skbio.util import get_data_path

from microprot.scripts.split_search import parse_pdb_match


class ParsersTests(TestCase):
    def setUp(self):
        self.file_fused = get_data_path('test_split_search/fused_fields.out')
        self.res_1 = {'Aligned_cols': 28,
                      'Cols': 28,
                      'E-value': 1.2e-07,
                      'Hit': ('3b85_A Phosphate starvation-inducible protein; '
                              'PHOH2, ATPase, PFAM: PF02562, ST genomics, PSI-'
                              '2, protein structure initiative; 2.35A {Coryneb'
                              'acterium glutamicum atcc 13032}'),
                      'Identities': 0.29,
                      'No': 1,
                      'P-value': 3.4e-12,
                      'Probab': 97.61,
                      'SS': 0.0,
                      'Score': 75.62,
                      'Similarity': 0.363,
                      'Sum_probs': 22.9,
                      'Template_Neff': 10.9,
                      'alignment': {
                          'Confidence': {'sequence':
                                         '4678999999999999999866665543'},
                          'Q Consensus': {'end': 250,
                                          'sequence':
                                          'pklifvqgaagtgktvllshlfyriaae',
                                          'start': 223,
                                          'totallen': 695},
                          'Q NZ_GG666849.1_': {'end': 250,
                                               'sequence':
                                               'PKLIFVQGAAGTGKTVLLSHLFYRIAAE',
                                               'start': 223,
                                               'totallen': 695},
                          'T 3b85_A': {'end': 49,
                                       'sequence':
                                       'NTIVFGLGPAGSGKTYLAMAKAVQALQS',
                                       'start': 22,
                                       'totallen': 208},
                          'T Consensus': {'end': 49,
                                          'sequence':
                                          '~~~~~i~g~~GtGKT~~~~~~~~~~~~~',
                                          'start': 22,
                                          'totallen': 208},
                          'T ss_dssp': {'sequence':
                                        'CSEEEEECCTTSSTTHHHHHHHHHHHHT'},
                          'T ss_pred': {'sequence':
                                        'CCeEEEECCCCCCHHHHHHHHHHHHHHh'},
                          'column score': {'sequence':
                                           '++.++++|+||||||+++..+..+++..'}}}

        self.res_362 = {
            'Aligned_cols': 28,
            'Cols': 28,
            'E-value': 0.022,
            'Hit': ('3cmu_A Protein RECA, recombinase A; homologous recombinat'
                    'ion, recombination/DNA complex; HET: DNA ADP; 4.20A {Esch'
                    'erichia coli}'),
            'Identities': 0.18,
            'No': 362,
            'P-value': 5.9e-07,
            'Probab': 92.19,
            'SS': 0.0,
            'Score': 62.61,
            'Similarity': 0.293,
            'Sum_probs': 24.1,
            'Template_Neff': 10.0,
            'alignment': {
                'Confidence': {'sequence': '5678999999999999999988777643'},
                'Q Consensus': {'end': 249,
                                'sequence': 'hpklifvqgaagtgktvllshlfyriaa',
                                'start': 222,
                                'totallen': 695},
                'Q NZ_GG666849.1_': {'end': 249,
                                     'sequence':
                                     'HPKLIFVQGAAGTGKTVLLSHLFYRIAA',
                                     'start': 222,
                                     'totallen': 695},
                'T 3cmu_A': {'end': 1107,
                             'sequence': 'MGRIVEIYGPESSGKTTLTLQVIAAAQR',
                             'start': 1080,
                             'totallen': 2050},
                'T Consensus': {'end': 1107,
                                'sequence': '~g~~ill~G~~G~GKT~la~~la~~~~~',
                                'start': 1080,
                                'totallen': 2050},
                'T ss_dssp': {'sequence': 'TTSEEEEECCTTSSHHHHHHHHHHHHHT'},
                'T ss_pred': {'sequence': 'CCeEEEEECCCCCCHHHHHHHHHHHHHH'},
                'column score': {'sequence': '.++++.+.|+||+|||+|+..+....+.'}}}

        self.res_454 = {
            'Aligned_cols': 25,
            'Cols': 25,
            'E-value': 0.052,
            'Hit': ('3j16_B RLI1P; ribosome recycling, translation, eukarya, '
                    'ribosome; HET: ATP; 7.20A {Saccharomyces cerevisiae} PDB'
                    ': 4crm_P*'),
            'Identities': 0.32,
            'No': 454,
            'P-value': 1.4e-06,
            'Probab': 90.39,
            'SS': 0.0,
            'Score': 54.77,
            'Similarity': 0.53,
            'Sum_probs': 20.9,
            'Template_Neff': 8.1,
            'alignment': {
                'Confidence': {'sequence': '4688999999999999998766554'},
                'Q Consensus': {'end': 248,
                                'sequence': 'klifvqgaagtgktvllshlfyria',
                                'start': 224,
                                'totallen': 695},
                'Q NZ_GG666849.1_': {'end': 248,
                                     'sequence': 'KLIFVQGAAGTGKTVLLSHLFYRIA',
                                     'start': 224,
                                     'totallen': 695},
                'T 3j16_B': {'end': 403,
                             'sequence': 'EILVMMGENGTGKTTLIKLLAGALK',
                             'start': 379,
                             'totallen': 608},
                'T Consensus': {'end': 403,
                                'sequence': 'ei~~i~G~nGsGKSTllk~l~G~~~',
                                'start': 379,
                                'totallen': 608},
                'T ss_dssp': {'sequence': 'CEEEEESCTTSSHHHHHHHHHTSSC'},
                'T ss_pred': {'sequence': 'cEEEEECCCCCCHHHHHHHHhcCCC'},
                'column score': {'sequence': '..+.+.|+.|+|||||++.+.-.+.'}}}

        self.res_455 = {
            'Aligned_cols': 26,
            'Cols': 26,
            'E-value': 0.052,
            'Hit': ('4f4c_A Multidrug resistance protein PGP-1; ABC transporte'
                    'r, ATPase, multi-drug transporter, exporter, A binding, h'
                    'ydrolase,protein transport; HET: NDG NAG BMA MAN 0SA; 3.4'
                    '0A {Caenorhabditis elegans}'),
            'Identities': 0.19,
            'No': 455,
            'P-value': 1.4e-06,
            'Probab': 90.38,
            'SS': 0.0,
            'Score': 56.06,
            'Similarity': 0.308,
            'Sum_probs': 21.2,
            'Template_Neff': 10.5,
            'alignment': {
                'Confidence': {'sequence': '46788999999999999988765543'},
                'Q Consensus': {'end': 249,
                                'sequence': 'klifvqgaagtgktvllshlfyriaa',
                                'start': 224,
                                'totallen': 695},
                'Q NZ_GG666849.1_': {'end': 249,
                                     'sequence': 'KLIFVQGAAGTGKTVLLSHLFYRIAA',
                                     'start': 224,
                                     'totallen': 695},
                'T 4f4c_A': {'end': 1131,
                             'sequence': 'QTLALVGPSGCGKSTVVALLERFYDT',
                             'start': 1106,
                             'totallen': 1321},
                'T Consensus': {'end': 1131,
                                'sequence': 'e~vaIvG~sGsGKSTLl~lL~g~~~~',
                                'start': 1106,
                                'totallen': 1321},
                'T ss_dssp': {'sequence': 'CEEEEECSTTSSTTSHHHHHTTSSCC'},
                'T ss_pred': {'sequence': 'CEEEEECCCCCCHHHHHHHHHhccCC'},
                'column score': {'sequence': '+.+.+.|+.|+|||||++.|.-....'}}}

        self.res_456 = {
            'Aligned_cols': 26,
            'Cols': 26,
            'E-value': 0.052,
            'Hit': ('3tqc_A Pantothenate kinase; biosynthesis of cofactors, pr'
                    'osthetic groups, carriers, TRAN; HET: ADP; 2.30A {Coxiell'
                    'a burnetii}'),
            'Identities': 0.23,
            'No': 456,
            'P-value': 1.4e-06,
            'Probab': 90.37,
            'SS': 0.0,
            'Score': 48.8,
            'Similarity': 0.271,
            'Sum_probs': 22.4,
            'Template_Neff': 9.3,
            'alignment': {
                'Confidence': {'sequence': '45788999999999999999988654'},
                'Q Consensus': {'end': 246,
                                'sequence': 'rhpklifvqgaagtgktvllshlfyr',
                                'start': 221,
                                'totallen': 695},
                'Q NZ_GG666849.1_': {'end': 246,
                                     'sequence': 'RHPKLIFVQGAAGTGKTVLLSHLFYR',
                                     'start': 221,
                                     'totallen': 695},
                'T 3tqc_A': {'end': 115,
                             'sequence': 'KVPYIIGIAGSVAVGKSTTSRVLKAL',
                             'start': 90,
                             'totallen': 321},
                'T Consensus': {'end': 115,
                                'sequence': '~~~~~i~i~G~~gsGKstl~~~l~~~',
                                'start': 90,
                                'totallen': 321},
                'T ss_dssp': {'sequence': 'CCCEEEEEECCTTSSHHHHHHHHHHH'},
                'T ss_pred': {'sequence': 'CCCeEEEEECCCCCCHHHHHHHHHHH'},
                'column score': {'sequence': '..+.+|.+.|++|+|||||++.|...'}}}

    def test_full(self):
        obs_full = parse_pdb_match(self.file_fused)
        exp_full = [self.res_1, self.res_362, self.res_454, self.res_455,
                    self.res_456]
        self.maxDiff = None
        self.assertEqual(obs_full, exp_full)


if __name__ == '__main__':
    main()
