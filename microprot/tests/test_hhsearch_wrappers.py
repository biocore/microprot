from unittest import TestCase, main
import os

from skbio.util import get_data_path

from microprot.scripts.pdb_search import mask_sequence


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

    def test_mask_sequence_filtering(self):
        s1 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_1-461'),
              ('MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALP'
               'NISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALIC'
               'RGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVL'
               'MAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKS'
               'MSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGIS'
               'NLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAER'
               'AMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGS'
               'SERSISVVVNNDDATTGVRVTHQMLFN'))
        s2 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_462-463'), 'TD')
        s3 = (('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;part'
               'ial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc'
               '_cont=0.531_464-815'),
              ('QVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQEEL'
               'AQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNKKANTSSMD'
               'YYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDE'
               'GMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEFNA'
               'EGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKN'
               'GENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLS'))
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
        exp = {'match': [s1, s3], 'non_match': [s2, s4]}
        self.assertEqual(obs, exp)

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
        seq = ('MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALP'
               'NISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALIC'
               'RGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVL'
               'MAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKS'
               'MSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGIS'
               'NLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAER'
               'AMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGS'
               'SERSISVVVNNDDATTGVRVTHQMLFN')
        header = ('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_2;'
                  'partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5'
                  '-10bp;gc_cont=0.531_1-461')
        header_nm = ('gi|556503834|ref|NC_000913.3|_2 # 337 # 2799 # 1 # ID=1_'
                     '2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spac'
                     'er=5-10bp;gc_cont=0.531_462-820')
        seq_nm = ('TDQVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLEN'
                  'WQEELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNK'
                  'KANTSSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSG'
                  'SLSYIFGKLDEGMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELA'
                  'DIEIEPVLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCR'
                  'VKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLS'
                  'WKLGV')

        exp = (header, seq)
        obs = mask_sequence(self.file_a, self.file_query,
                            min_fragment_length=450)
        self.assertEqual(obs['match'][0], exp)

        filename = '/tmp/test.mfa'
        mask_sequence(self.file_a, self.file_query, filename,
                      min_fragment_length=450)

        f = open(filename+'.match', 'r')
        obs = f.readlines()
        f.close()
        os.remove(filename+'.match')
        self.assertIn(seq+"\n", obs)
        self.assertIn(">"+header+"\n", obs)

        f = open(filename+'.non_match', 'r')
        obs = f.readlines()
        f.close()
        os.remove(filename+'.non_match')
        self.assertIn(seq_nm+"\n", obs)
        self.assertIn(">"+header_nm+"\n", obs)

        with self.assertRaises(IOError):
            mask_sequence(self.file_a, self.file_query, '/dev')


if __name__ == '__main__':
    main()
