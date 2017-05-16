from unittest import TestCase, main
import os
import sys

from filecmp import cmp
from shutil import rmtree
from io import StringIO
from contextlib import contextmanager

from skbio.util import get_data_path

from microprot.scripts.split_search import mask_sequence, pretty_output


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class SplitSeq(TestCase):
    def setUp(self):
        self.file_a = get_data_path('test_split_search/NC_000913.3_2.out')
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
        self.file_query = get_data_path(
                          'test_split_search/NC_000913.3_2.fasta')

    def test_mask_sequence_filtering(self):
        s1 = (('gi|556503834|ref|NC_000913.3|_2_1-461 # 337 # 2799 # 1 # ID=1_'
               '2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer='
               '5-10bp;gc_cont=0.531'),
              ('MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALP'
               'NISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALIC'
               'RGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVL'
               'MAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKS'
               'MSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGIS'
               'NLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAER'
               'AMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGS'
               'SERSISVVVNNDDATTGVRVTHQMLFN'))
        s2 = (('gi|556503834|ref|NC_000913.3|_2_462-463 # 337 # 2799 # 1 # ID='
               '1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer='
               '5-10bp;gc_cont=0.531'), 'TD')
        s3 = (('gi|556503834|ref|NC_000913.3|_2_464-815 # 337 # 2799 # 1 # ID='
               '1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer='
               '5-10bp;gc_cont=0.531'),
              ('QVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQEEL'
               'AQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNKKANTSSMD'
               'YYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDE'
               'GMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEFNA'
               'EGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKN'
               'GENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLS'))
        s4 = (('gi|556503834|ref|NC_000913.3|_2_816-820 # 337 # 2799 # 1 # ID='
               '1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer='
               '5-10bp;gc_cont=0.531', 'WKLGV'))
        s5 = (('gi|556503834|ref|NC_000913.3|_2_462-820 # 337 # 2799 # 1 # ID'
               '=1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer='
               '5-10bp;gc_cont=0.531'),
              ('TDQVIEVFVIGVGGVGGALLEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQE'
               'ELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAVADQYADFLREGFHVVTPNKKANTSS'
               'MDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKL'
               'DEGMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEF'
               'NAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKV'
               'KNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLSWKLGV'))
        s6 = (('gi|556503834|ref|NC_000913.3|_2_1-820 # 337 # 2799 # 1 # ID='
               '1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer='
               '5-10bp;gc_cont=0.531'), self.query)

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
        self.assertEqual(obs, {'match': [s1, s3], 'non_match': []})

        obs = mask_sequence(self.file_a, self.file_query, min_prob=99.0,
                            max_evalue=4.90e-41, max_pvalue=0.00011,
                            min_fragment_length=4)
        self.assertEqual(obs, {'match': [s1, s3], 'non_match': [s4]})

    def test_mask_sequence_information(self):
        seq = ('MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALP'
               'NISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALIC'
               'RGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVL'
               'MAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKS'
               'MSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGIS'
               'NLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAER'
               'AMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGS'
               'SERSISVVVNNDDATTGVRVTHQMLFN')
        header = ('gi|556503834|ref|NC_000913.3|_2_1-461 # 337 # 2799 # 1 # ID'
                  '=1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spac'
                  'er=5-10bp;gc_cont=0.531')

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
        self.assertFalse(obs)

        with self.assertRaises(IOError):
            mask_sequence(self.file_a, self.file_query, '/dev')

    def test_pretty_output(self):
        pretty_fp = get_data_path('test_split_search/NC_000913.3_2.pretty')
        with open(pretty_fp, 'r') as f:
            pretty = f.read()
        mask_obs = mask_sequence(self.file_a, self.file_query,
                                 min_fragment_length=5)
        with captured_output() as (out, err):
            pretty_output(mask_obs)
        output = out.getvalue()
        self.assertEqual(output, pretty)


if __name__ == '__main__':
    main()
