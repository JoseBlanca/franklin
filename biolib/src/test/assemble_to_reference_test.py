'''
Created on 15/09/2009

@author: peio
'''
import unittest
from StringIO import StringIO
from tempfile import NamedTemporaryFile
from biolib.libassemble import check_cfile, get_files_to_load_in_bank
from biolib.biolib_utils import NamedTemporaryDir
import biolib
import os.path

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class AssembleToReferenceTest(unittest.TestCase):
    'This class tests functions in assemble to reference'

    def test_check_cfile(self):
        'It test check_cfile'
        cfhand = StringIO( """{'work_dir'  : 'working_dir',
                 'reference' : [{'name'      : 'reference',
                                  'seq_fpath' : 'reference.fasta',
                                  'format':'fasta'}],
                 'reads'      : [{'name'       :'solexa1',
                                 'seq_fpath'  : 'name1',
                                 'qual_fpath' : 'name1',
                                 'format'     : 'format1',
                                 'seq_type'   : 'solexa',
                                 'aligner'    : 'nucmer',
                                 'aligner_parameters' : '-l 20 -c 20'}]}""")
        configuration = check_cfile(cfhand)
        assert configuration['reads'][0]['name'] == 'solexa1'
        assert configuration['reference'][0]['qual_fpath'] == None

        ## Without requiered fields
        cfhand = StringIO( """{}""")
        try:
            configuration = check_cfile(cfhand)
            self.fail()
            #pylint: disable-msg=W0704
        except RuntimeError:
            pass
        ###
        ref_fhand  = open(os.path.join(DATA_DIR, 'seq.fasta'), 'r')
        read_fhand = open(os.path.join(DATA_DIR, 'seq.fasta'), 'r')
        dict_to_eval =  """{'work_dir'  : 'working_dir',
                             'reference' : [{'name'      : 'reference',
                             'seq_fpath' : '""" + ref_fhand.name+ """'}],
                              'reads'     : [{'name'       :'sanger',
                              'seq_fpath'  : '""" + read_fhand.name +"""'}]}"""
        cfhand = StringIO(dict_to_eval)
        configuration = check_cfile(cfhand)
        assert configuration['reads'][0]['aligner'] == 'nucmer'
        assert configuration['reads'][0]['format']  == 'fasta'
        assert configuration['reference'][0]['format']  == 'fasta'

    @staticmethod
    def test_get_files_to_load_in_bank():
        'It tests get_files_to_load_in_bank'
        temp_dir = NamedTemporaryDir()
        ref_fhand  = open(os.path.join(DATA_DIR, 'seq.fasta'), 'r')
        read_fhand = NamedTemporaryFile(suffix='.fastq')
        read_fhand.write('''@SRR019165.1 :5:1:898:110 length=25
TGTAAGGGAGCAGCGGAGTGGGCAA
+SRR019165.1 :5:1:898:110 length=25
IIIIIIIIIIIII+*IIIII.II5I
@SRR019165.2 :5:1:902:102 length=25
TAGTACCTCATTCTGGTTTTAATTG
+SRR019165.2 :5:1:902:102 length=25
IIIIIIHIIIIIIIIIIIIIIIIII
@SRR019165.3 :5:1:941:116 length=25
AGGTCGAGGCATGAGGCTGGAGCAC
+SRR019165.3 :5:1:941:116 length=25
IIIIIIIIIIIIII*IAIIIII.II''')
        read_fhand.flush()
        dict_to_eval =  """{'work_dir'  : '"""+temp_dir.get_name()+"""',
                            'reference' : [{'name'      : 'reference',
                            'seq_fpath' : '""" + ref_fhand.name+ """'}],
                            'reads'     : [{'name'       :'sanger',
                            'seq_fpath'  : '""" + read_fhand.name +"""',
                            'format':'fastq'}]}"""
        cfhand = StringIO(dict_to_eval)
        configuration = check_cfile(cfhand)
        files = get_files_to_load_in_bank(configuration)
        assert len(files) == 2
        assert files[0] == os.path.join(DATA_DIR, 'seq.fasta')



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()