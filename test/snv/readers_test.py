'''
Created on 2011 aza 21

@author: peio
'''
import unittest
from os.path import join
from franklin.utils.misc_utils import TEST_DATA_DIR
from franklin.snv.readers import VcfParser

class TestVcfParser(unittest.TestCase):


    def test_vcfparser(self):
        'It test the vcf arser'
        vcf_path = join(TEST_DATA_DIR, 'contigs.vcf')
        vcf = VcfParser(vcf_path)
        assert vcf.version == '3.3'

        header = vcf.header
        filter_info = 'SNV is not a CAP detectable by the enzymes: cheap ones'
        assert header['FILTER']['CEF'] == filter_info

        vcf_ = vcf.get_snv(('CUTC000002', '79'))
        assert vcf_['CHROM'] == 'CUTC000002'
        assert vcf_['POS'] == '79'

        vcfs = list(vcf.vcfs)
        assert len(vcfs) == 15
        assert  vcfs[4]['samples'] == {'MU16_454_MU16': {},
                                       'UPV196_454_UPV196': {'C': 2, 'T': 16}}
        assert vcfs[14]['samples'] == {'MU16_454_MU16': {'A': 2},
                                       'UPV196_454_UPV196': {'G': 2}}

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_vcfparser']
    unittest.main()
