from unittest import TestCase, main

from skbio.util import get_data_path

from microprot.scripts.pdb_search import is_overlapping, \
                                         _parse_hit_summary_line, \
                                         _parse_hit_block, parse_pdb_match, \
                                         select_hits, report_hits, \
                                         report_uncovered_subsequences


class ParsersTests(TestCase):
    def setUp(self):
        self.line_a = ('  3 2cdq_A Aspartokinase; aspartat 100.0 4.1E-58 '
                       '1.1E-62  453.8   0.0  458    1-473    28-496 (510)\n')
        self.true_a = {'No': 3,
                       'Hit': '2cdq_A Aspartokinase; aspartat',
                       'Prob': 100.0,
                       'E-value': 4.1e-58,
                       'P-value': 1.1e-62,
                       'Score': 453.8,
                       'SS': 0.0,
                       'Cols': 458,
                       'Query HMM': '1-473',
                       'Template HMM': '28-496',
                       '#match states': 510,
                       }
        self.line_b = ('    36 1gtm_A Glutamate dehydrogenase  98.8 2.9E-12 '
                       '7.9E-17  126.1   0.0  201  451-679   193-407 (419)\n')
        self.true_b = {'No': 36,
                       'Hit': '1gtm_A Glutamate dehydrogenase',
                       'Prob': 98.8,
                       'E-value': 2.9E-12,
                       'P-value': 7.9E-17,
                       'Score': 126.1,
                       'SS': 0.0,
                       'Cols': 201,
                       'Query HMM': '451-679',
                       'Template HMM': '193-407',
                       '#match states': 419,
                       }
        self.line_c = ('  5 3hdi_A Processing protease; CA  16.4      57  '
                       '0.0016   23.0   0.0   81   20-100   124-206 (421)')
        self.true_c = {'No': 5,
                       'Hit': '3hdi_A Processing protease; CA',
                       'Prob': 16.4,
                       'E-value': 57,
                       'P-value': 0.0016,
                       'Score': 23.0,
                       'SS': 0.0,
                       'Cols': 81,
                       'Query HMM': '20-100',
                       'Template HMM': '124-206',
                       '#match states': 421,
                       }
        self.line_d = ('  9 1hr6_B Beta-MPP, mitochondrial  12.5      84  '
                       '0.0023   22.0   0.0   81   20-100   129-212 (443)\n')
        self.true_d = {'No': 9,
                       'Hit': '1hr6_B Beta-MPP, mitochondrial',
                       'Prob': 12.5,
                       'E-value': 84,
                       'P-value': 0.0023,
                       'Score': 22.0,
                       'SS': 0.0,
                       'Cols': 81,
                       'Query HMM': '20-100',
                       'Template HMM': '129-212',
                       '#match states': 443,
                       }
        self.block_a = ("No 1\n"
                        ">2xqo_A Cellulosome enzyme, dockerin type I; hydrolas"
                        "e; HET: MSE CTR; 1.40A {Clostridium thermocellum}\n"
                        "Probab=31.45  E-value=18  Score=28.87  Aligned_cols=3"
                        "8  Identities=29%  Similarity=0.348  Sum_probs=30.9  "
                        "Template_Neff=1.400\n"
                        "\n"
                        "Q T0810-D1         71 ENEIRDENNRALVSGCKEELLAQVDEHFDEQ"
                        "RESMSVP  108 (113)\n"
                        "Q Consensus        71 eneirdennralvsgckeellaqvdehfdeq"
                        "resmsvp  108 (113)\n"
                        "                      ..|||.-|...-|+..|-.|..|||||+|-."
                        "+.....|\n"
                        "T Consensus        43 ~~Ei~A~~~~m~VqEVK~~l~KeiDEHWdlI"
                        "~~~~G~~   80 (243)\n"
                        "T 2xqo_A           43 KHELIARASSLKVSEVKAIIKKQVDEHWDVI"
                        "RDVCGFK   80 (243)\n"
                        "T ss_dssp             HHHHHHHHHHCCHHHHHHHHHHHHHHTHHHH"
                        "HHHHCCS\n"
                        "T ss_pred             HHHHHHHhccccHHHHHHHHHHHHHHHHHHH"
                        "HHHhCCC\n"
                        "Confidence            3578888888899999999999999999976"
                        "6544433\n")
        self.true_block_a = {
            'No': 1,
            'Hit': ('2xqo_A Cellulosome enzyme, dockerin type I; hydrolase; H'
                    'ET: MSE CTR; 1.40A {Clostridium thermocellum}'),
            'Probab': 31.45,
            'E-value': 18.0,
            'Score': 28.87,
            'Aligned_cols': 38,
            'Identities': 0.29,
            'Similarity': 0.348,
            'Sum_probs': 30.9,
            'Template_Neff': 1.4,
            'alignment': {
                'Q T0810-D1': {'start': 71, 'end': 108, 'totallen': 113,
                               'sequence': ('ENEIRDENNRALVSGCKEELLAQVDEHFDEQR'
                                            'ESMSVP')},
                'Q Consensus': {'start': 71, 'end': 108, 'totallen': 113,
                                'sequence': ('eneirdennralvsgckeellaqvdehfdeq'
                                             'resmsvp')},
                'column score': {'sequence': ('..|||.-|...-|+..|-.|..|||||+|-'
                                              '.+.....|')},
                'T Consensus': {'start': 43, 'end': 80, 'totallen': 243,
                                'sequence': ('~~Ei~A~~~~m~VqEVK~~l~KeiDEHWdlI'
                                             '~~~~G~~')},
                'T 2xqo_A': {'start': 43, 'end': 80, 'totallen': 243,
                             'sequence': ('KHELIARASSLKVSEVKAIIKKQVDEHWDVIRDV'
                                          'CGFK')},
                'T ss_dssp': {'sequence': ('HHHHHHHHHHCCHHHHHHHHHHHHHHTHHHHHH'
                                           'HHCCS')},
                'T ss_pred': {'sequence': ('HHHHHHHhccccHHHHHHHHHHHHHHHHHHHHH'
                                           'HhCCC')},
                'Confidence': {'sequence': ('35788888888999999999999999999766'
                                            '544433')}
            }}
        self.block_b = (
            "No 10\n"
            ">4xea_A Peptidase M16 domain protein; metallopeptidase, structur"
            "al genomics, PSI-biology, midwest for structural genomics, MCSG;"
            " 1.95A {Alicyclobacillus acidocaldarius subsp}\n"
            "Probab=12.10  E-value=89  Score=21.62  Aligned_cols=80  Identiti"
            "es=13%  Similarity=0.110  Sum_probs=43.1  Template_Neff=10.600\n"
            "\n"
            "Q T0810-D1         21 VEWENSEANPEALFANWRHEFMVDSSKRESMKTELCKELQAL"
            "PAQDLTLFENEIRDENNR-ALVSGC--KEELLAQVDEH   97 (113)\n"
            "Q Consensus        21 vewenseanpealfanwrhefmvdsskresmktelckelqal"
            "paqdltlfeneirdennr-alvsgc--keellaqvdeh   97 (113)\n"
            "                      -|+.....+|+.++..--++.+-..+....-..--.+.+..+"
            "...|+.-|-+..-..+|- ..+.|.  .+++...+.++\n"
            "T Consensus       135 ~ei~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~G~~~~l~~i"
            "t~~~l~~~~~~~y~~~~~~lvi~Gd~~~~~~~~l~~~~  214 (423)\n"
            "T 4xea_A          135 QEIHMVNDHPDRRAYMELLRAMYHEHPVRIDIAGTVESVRAI"
            "TKEQLLLCYDTFYHPSNMVLVIAGGFDADEIAHVIEEN  214 (423)\n"
            "T ss_dssp             HHHHHHHSCHHHHHHHHHHHHHCSSCGGGSCTTCCHHHHHHC"
            "CHHHHHHHHHHHCSGGGEEEEEEESSCHHHHHHHHHHH\n"
            "T ss_pred             HHHHHhcCChHHHHHHHHHHHhcccCCCCCCCCCCHHHHhhC"
            "CHHHHHHHHHHhCCccceEEEEECCCCHHHHHHHHHHh\n"
            "Confidence            345455556665543322333332221111111234567888"
            "888888887776655444 456666  45677777776\n"
            "\n"
            "\n"
            "Q T0810-D1         98 FDE  100 (113)\n"
            "Q Consensus        98 fde  100 (113)\n"
            "                      |..\n"
            "T Consensus       215 ~~~  217 (423)\n"
            "T 4xea_A          215 QAK  217 (423)\n"
            "T ss_dssp             HHT\n"
            "T ss_pred             hcc\n"
            "Confidence            653\n"
        )
        self.true_block_b = {
            'Score': 21.62,
            'E-value': 89.0,
            'Template_Neff': 10.6,
            'Probab': 12.1,
            'Sum_probs': 43.1,
            'Similarity': 0.11,
            'No': 10,
            'alignment': {
                'column score': {'sequence': ('-|+.....+|+.++..--++.+-..+....-'
                                              '..--.+.+..+...|+.-|-+..-..+|- .'
                                              '.+.|.  .+++...+.++|..')},
                'T Consensus': {'start': 135, 'end': 217, 'totallen': 423,
                                'sequence': ('~ei~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                                             '~G~~~~l~~it~~~l~~~~~~~y~~~~~~lvi'
                                             '~Gd~~~~~~~~l~~~~~~~')},
                'Q T0810-D1': {'start': 21, 'end': 100, 'totallen': 113,
                               'sequence': ('VEWENSEANPEALFANWRHEFMVDSSKRESMKT'
                                            'ELCKELQALPAQDLTLFENEIRDENNR-ALVSG'
                                            'C--KEELLAQVDEHFDE')},
                'T ss_dssp': {'sequence': ('HHHHHHHSCHHHHHHHHHHHHHCSSCGGGSCTTC'
                                           'CHHHHHHCCHHHHHHHHHHHCSGGGEEEEEEESS'
                                           'CHHHHHHHHHHHHHT')},
                'Confidence': {'sequence': ('345455556665543322333332221111111'
                                            '234567888888888887776655444 45666'
                                            '6  45677777776653')},
                'T 4xea_A': {'start': 135, 'end': 217, 'totallen': 423,
                             'sequence': ('QEIHMVNDHPDRRAYMELLRAMYHEHPVRIDIAGT'
                                          'VESVRAITKEQLLLCYDTFYHPSNMVLVIAGGFDA'
                                          'DEIAHVIEENQAK')},
                'T ss_pred': {'sequence': ('HHHHHhcCChHHHHHHHHHHHhcccCCCCCCCCC'
                                           'CHHHHhhCCHHHHHHHHHHhCCccceEEEEECCC'
                                           'CHHHHHHHHHHhhcc')},
                'Q Consensus': {'start': 21, 'end': 100, 'totallen': 113,
                                'sequence': ('vewenseanpealfanwrhefmvdsskresmk'
                                             'telckelqalpaqdltlfeneirdennr-alv'
                                             'sgc--keellaqvdehfde')}},
            'Aligned_cols': 80,
            'Hit': ('4xea_A Peptidase M16 domain protein; metallopeptidase, st'
                    'ructural genomics, PSI-biology, midwest for structural ge'
                    'nomics, MCSG; 1.95A {Alicyclobacillus acidocaldarius subs'
                    'p}'),
            'Identities': 0.13}
        self.file_a = get_data_path('test_pdb_search/NC_000913.3_2.out')
        self.file_b = get_data_path('test_pdb_search/T0810-D1.fasta.out')
        self.true_a_hits_report = [{
            'end': 461,
            'pdb_id': '2j0w_A',
            'start': 1,
            'covered_sequence':
                ('MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDA'
                 'LPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINA'
                 'ALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIP'
                 'ADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQV'
                 'PDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRD'
                 'EDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISF'
                 'CVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAAL'
                 'ARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFN')}, {
            'end': 815,
            'pdb_id': '1ebf_A',
            'start': 464,
            'covered_sequence':
                ('QVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQE'
                 'ELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNKKANT'
                 'SSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYI'
                 'FGKLDEGMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEP'
                 'VLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVD'
                 'GNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLS')}]
        self.query = (
            'MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNIS'
            'DAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMS'
            'IAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNE'
            'KGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYF'
            'GAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGM'
            'KGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPL'
            'AVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVT'
            'HQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLEN'
            'WQEELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNKKANTSS'
            'MDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDEG'
            'MSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEFNAEGDV'
            'AAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKNGENALAF'
            'YSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLSWKLGV')
        self.true_subseqs = [{'sequence': 'TD', 'start': 462, 'end': 463},
                             {'sequence': 'WKLGV', 'start': 816, 'end': 820}]

    def test_is_overlapping(self):
        self.assertTrue(is_overlapping((20, 100), (30, 150)))
        self.assertTrue(is_overlapping((20, 100), (2, 70)))
        self.assertTrue(is_overlapping((20, 100), (30, 70)))
        self.assertTrue(is_overlapping((20, 100), (5, 150)))
        self.assertFalse(is_overlapping((20, 100), (150, 200)))
        self.assertFalse(is_overlapping((20, 100), (1, 8)))

        with self.assertRaises(ValueError):
            is_overlapping((100, 20), (1, 8))
        with self.assertRaises(ValueError):
            is_overlapping((10, 20), (-9, -20))

    def test__parse_hit_summary_line(self):
        self.assertEqual(self.true_a, _parse_hit_summary_line(self.line_a))
        self.assertEqual(self.true_b, _parse_hit_summary_line(self.line_b))
        self.assertEqual(self.true_c, _parse_hit_summary_line(self.line_c))
        self.assertEqual(self.true_d, _parse_hit_summary_line(self.line_d))
        with self.assertRaises(ValueError):
            _parse_hit_summary_line('no valid line')

    def test__parse_hit_block(self):
        self.assertEqual(self.true_block_a, _parse_hit_block(self.block_a))
        self.assertEqual(self.true_block_b, _parse_hit_block(self.block_b))

    def test_parse_pdb_match(self):
        res = parse_pdb_match(self.file_b)
        self.assertEqual(len(res), 10)

        r0 = res[0]
        del r0['Cols']
        del r0['P-value']
        del r0['SS']
        self.assertEqual(self.true_block_a, r0)

        r9 = res[9]
        del r9['Cols']
        del r9['P-value']
        del r9['SS']
        self.assertEqual(self.true_block_b, r9)

    def test_select_hits(self):
        hits = parse_pdb_match(self.file_b)
        self.assertEqual(len(select_hits(hits, e_value_threshold=100)), 2)
        self.assertEqual(len(select_hits(hits)), 0)

        hits = parse_pdb_match(self.file_a)
        goodhits = select_hits(hits)
        self.assertEqual(len(goodhits), 2)
        self.assertEqual(goodhits[0]['No'], 1)
        self.assertEqual(goodhits[1]['No'], 5)

    def test_report_hits(self):
        hits = select_hits(parse_pdb_match(self.file_a))
        rep = report_hits(hits)
        self.assertEqual(self.true_a_hits_report, rep)

    def test_report_uncovered_subsequences(self):
        hits = select_hits(parse_pdb_match(self.file_a))
        subseqs = report_uncovered_subsequences(hits, self.query, 0)
        self.assertEqual(subseqs, self.true_subseqs)

if __name__ == '__main__':
    main()
