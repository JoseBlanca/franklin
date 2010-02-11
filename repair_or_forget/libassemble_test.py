'''
Created on 15/09/2009

@author: peio
'''
import unittest
from franklin.libassemble import (check_and_fix_config,
                                prepare_files_to_load_in_bank)
from franklin.utils.misc_utils import NamedTemporaryDir
import franklin
import os.path

DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'data')

class AssembleToReferenceTest(unittest.TestCase):
    'This class tests functions in assemble to reference'

    def test_check_and_fix_config(self):
        'It test the config checker'
        temp_dir = NamedTemporaryDir()
        config = {'work_dir'  : temp_dir.get_name(),
                 'reference' : [{'name'      : 'reference',
                                  'seq_fpath' : 'reference.fasta',
                                  'format':'fasta'}],
                 'reads'      : [{'name'       :'solexa1',
                                 'seq_fpath'  : 'name1',
                                 'qual_fpath' : 'name1',
                                 'format'     : 'format1',
                                 'seq_type'   : 'solexa',
                                 'aligner'    : 'nucmer',
                                 'aligner_parameters' : '-l 20 -c 20'}]}
        check_and_fix_config(config)
        assert config['reads'][0]['name'] == 'solexa1'
        assert config['reference'][0]['qual_fpath'] == None

        ## Without requiered fields
        cfhand = {}
        try:
            check_and_fix_config(cfhand)
            self.fail()
            #pylint: disable-msg=W0704
        except ValueError:
            pass
        ### Only with required fields
        ref_fhand  = open(os.path.join(DATA_DIR, 'seq.fasta'), 'r')
        read_fhand = open(os.path.join(DATA_DIR, 'seq.fasta'), 'r')
        config =  {'work_dir'  : temp_dir.get_name(),
                   'reference' : [{'name'     : 'reference',
                                   'seq_fpath': ref_fhand.name}],
                   'reads'     : [{'name'     :'sanger',
                                   'seq_fpath': read_fhand.name}]}
        check_and_fix_config(config)
        assert config['reads'][0]['aligner'] == 'nucmer'
        assert config['reads'][0]['format']  == 'fasta'
        assert config['reference'][0]['format']  == 'fasta'

    @staticmethod
    def test_prepare_files_to_load_in_bank():
        'It tests get_files_to_load_in_bank'
        temp_dir = NamedTemporaryDir()
        ref_fhand  = open(os.path.join(DATA_DIR, 'seq.fasta'), 'r')
        read_fhand = open(os.path.join(DATA_DIR, 'solexa.fastq'), 'r')
        config =  {'work_dir' : temp_dir.get_name(),
                   'reference': [{'name'     : 'reference',
                                  'seq_fpath': ref_fhand.name}],
                    'reads'   : [{'name'     : 'sanger',
                                  'seq_fpath': read_fhand.name,
                                  'format'   :'fastq'}]}
        check_and_fix_config(config)
        files = prepare_files_to_load_in_bank(config)
        assert len(files) == 2
        assert files[0] == '%s/seq.fasta.seq' % temp_dir.get_name()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()