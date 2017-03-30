from unittest import TestCase, main
import os

from skbio.util import get_data_path

from microprot.scripts.pdb_search import mask_sequence, frag_size


class SplitSeq(TestCase):
    def setUp(self):
        self.file_a = get_data_path('test_pdb_search/NC_000913.3_2.out')
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
        self.file_query = get_data_path('test_pdb_search/NC_000913.3_2.fasta')

        self.hit = {
            'No': 225,
            'Hit': ('3ofg_A BOCA/MESD chaperone for YWTD beta-propeller-EGF P;'
                    ' molecular chaperone, protein folding, YWTD propeller, LD'
                    'LR; HET: MSE; 1.37A {Caenorhabditis elegans}'),
            'Similarity': 0.24,
            'alignment': {
                'Q gi|556503834|r': {
                    'start': 483, 'totallen': 820, 'end': 525,
                    'sequence': 'QLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQEEL'},
                'T Consensus': {
                    'start': 35, 'totallen': 95, 'end': 77,
                    'sequence': 'ia~~Wq~~L~n~~I~v~~y~vd~~r~if~~~dG~~a~e~k~FL'},
                'Q Consensus': {
                    'start': 483, 'totallen': 820, 'end': 525,
                    'sequence': 'qlkrqqswlknkhidlrvcgvanskalltnvhglnlenwqeel'},
                'Confidence': {
                    'sequence': '3456788899999999999999999999988887665444433'},
                'column score': {
                    'sequence': '--++-|+-|.|.||+.++.+|..++++.+--+|-......+=|'},
                'T 3ofg_A': {
                    'start': 35, 'totallen': 95, 'end': 77, 'sequence':
                    'WTQIWQSQLYNNHVDLQVFVIDDNRAIFMFKNGEQAFEAKKFL'},
                'T ss_dssp': {
                    'sequence': 'HHHHHHHHHHHTTCCEEEEEEETTEEEEEESSGGGHHHHHHHH'},
                'T ss_pred': {
                    'sequence': 'HHHHHHHHHHhCCeEEEEEEEcCCEEEEEeCchhhHHHHHHHH'}
            },
            'Aligned_cols': 43,
            'SS': 0.0,
            'Score': 28.65,
            'Probab': 20.11,
            'Cols': 43,
            'P-value': 0.0011,
            'Identities': 0.26,
            'E-value': 41.0,
            'Sum_probs': 34.5,
            'Template_Neff': 5.6}

        self.minhit = {
            'alignment': {
                'Q gi|556503834|r': {
                    'sequence': 'QLKRQQSWLKNK--ALL---T-NVHGLNLENWQEEL'}}
            }

    def test_mask_sequence_filtering(self):
        s1 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_1-461'),
              ('mrvlkfggtsvanaerflrvadilesnarqgqvatvlsapakitnhlvamiektisgqdalp'
               'nisdaerifaelltglaaaqpgfplaqlktfvdqefaqikhvlhgisllgqcpdsinaalic'
               'rgekmsiaimagvlearghnvtvidpvekllavghylestvdiaestrriaasripadhmvl'
               'magftagnekgelvvlgrngsdysaavlaaclradcceiwtdvdgvytcdprqvpdarllks'
               'msyqeamelsyfgakvlhprtitpiaqfqipclikntgnpqapgtligasrdedelpvkgis'
               'nlnnmamfsvsgpgmkgmvgmaarvfaamsrarisvvlitqssseysisfcvpqsdcvraer'
               'amqeefylelkeglleplavterlaiisvvgdgmrtlrgisakffaalaraninivaiaqgs'
               'sersisvvvnnddattgvrvthqmlfn'))
        s2 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_462-463'), 'TD')
        s3 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_464-815'),
              ('qvievfvigvggvggalleqlkrqqswlknkhidlrvcgvanskalltnvhglnlenwqeel'
               'aqakepfnlgrlirlvkeyhllnpvivdctssqavadqyadflregfhvvtpnkkantssmd'
               'yyhqlryaaeksrrkflydtnvgaglpvienlqnllnagdelmkfsgilsgslsyifgklde'
               'gmsfseattlaremgytepdprddlsgmdvarkllilaretgreleladieiepvlpaefna'
               'egdvaafmanlsqlddlfaarvakardegkvlryvgnidedgvcrvkiaevdgndplfkvkn'
               'genalafyshyyqplplvlrgygagndvtaagvfadllrtls'))
        s4 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_816-820', 'WKLGV'))
        s5 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_462-820'),
              ('TDQVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQE'
               'ELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNKKANTSS'
               'MDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKL'
               'DEGMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEF'
               'NAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKV'
               'KNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLSWKLGV'))
        s6 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_1-820'), self.query)

        obs = mask_sequence(self.file_a, self.file_query, None)
        self.assertEqual(obs, {'match': [s1, s3], 'non_match': [s2, s4]})

        obs = mask_sequence(self.file_a, self.file_query, min_prob=100.0)
        self.assertEqual(obs, {'match': [s1, s3], 'non_match': [s2, s4]})

        obs = mask_sequence(self.file_a, self.file_query, max_evalue=4.1e-58)
        self.assertEqual(obs, {'match': [s1], 'non_match': [s5]})

        obs = mask_sequence(self.file_a, self.file_query, max_pvalue=1e-58)
        self.assertEqual(obs, {'match': [s1], 'non_match': [s5]})

        obs = mask_sequence(self.file_a, self.file_query,
                            min_fragment_length=500)
        self.assertEqual(obs, {'match': [], 'non_match': [s6]})

        obs = mask_sequence(self.file_a, self.file_query, min_prob=99.0,
                            max_evalue=4.90e-41, max_pvalue=0.00011,
                            min_fragment_length=200)
        self.assertEqual(obs, {'match': [s1, s3], 'non_match': [s2, s4]})

    def test_mask_sequence_information(self):
        seq = ('mrvlkfggtsvanaerflrvadilesnarqgqvatvlsapakitnhlvamiektisgqdalp'
               'nisdaerifaelltglaaaqpgfplaqlktfvdqefaqikhvlhgisllgqcpdsinaalic'
               'rgekmsiaimagvlearghnvtvidpvekllavghylestvdiaestrriaasripadhmvl'
               'magftagnekgelvvlgrngsdysaavlaaclradcceiwtdvdgvytcdprqvpdarllks'
               'msyqeamelsyfgakvlhprtitpiaqfqipclikntgnpqapgtligasrdedelpvkgis'
               'nlnnmamfsvsgpgmkgmvgmaarvfaamsrarisvvlitqssseysisfcvpqsdcvraer'
               'amqeefylelkeglleplavterlaiisvvgdgmrtlrgisakffaalaraninivaiaqgs'
               'sersisvvvnnddattgvrvthqmlfn')
        header = ('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;'
                  'partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5'
                  '-10bp;gc_cont=0.531_1-461')
        exp = (header, seq)
        obs = mask_sequence(self.file_a, self.file_query,
                            min_fragment_length=450)
        self.assertEqual(obs['match'][0], exp)

        filename = '/tmp/test.mfa'
        mask_sequence(self.file_a, self.file_query, filename,
                      min_fragment_length=450)

        f = open(filename, 'r')
        obs = f.readlines()
        f.close()
        os.remove(filename)
        self.assertIn(seq+"\n", obs)
        self.assertIn(">"+header+"\n", obs)

    def test_frag_size(self):
        obs = frag_size(self.minhit)
        self.assertEqual(obs, 30)

        obs = frag_size(self.hit)
        self.assertEqual(obs, 43)

        self.assertRaisesRegex(KeyError,
                               'alignment',
                               frag_size,
                               {})

        self.assertRaisesRegex(IndexError,
                               'list index out of range',
                               frag_size,
                               {'alignment': {}})

        self.assertRaisesRegex(KeyError,
                               'sequence',
                               frag_size,
                               {'alignment': {'Q consensus': {}}})


if __name__ == '__main__':
    main()
