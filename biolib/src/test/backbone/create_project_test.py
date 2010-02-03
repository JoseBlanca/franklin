'''
Created on 26/01/2010

@author: jose
'''
import unittest, os.path

from biolib.utils.misc_utils import NamedTemporaryDir
from biolib.backbone.create_project import create_project
from biolib.backbone.analysis import do_analysis, BACKBONE_DIRECTORIES
from configobj import ConfigObj

class TestBackbone(unittest.TestCase):
    'It tests the backbone'

    @staticmethod
    def test_create_project():
        'We can create a project'
        test_dir = NamedTemporaryDir()
        settings_path = create_project(directory=test_dir.name,
                                       name='backbone')

        assert settings_path == os.path.join(test_dir.name,
                                'backbone', BACKBONE_DIRECTORIES['config_file'])
        settings = ConfigObj(settings_path)
        assert settings['General_settings']['project_name'] == 'backbone'
        project_path = os.path.join(test_dir.name, 'backbone')
        assert settings['General_settings']['project_path'] == project_path
        test_dir.close()

    @staticmethod
    def test_cleaning_analysis():
        'We can clean the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name)
        project_dir = os.path.join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = os.path.join(project_dir, 'reads')
        original_reads_dir = os.path.join(reads_dir, 'original')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)

        #print original_reads_dir
        fpath_454 = os.path.join(original_reads_dir, 'pt_454.sfastq')
        reads_454 = '@seq1\nGCTAGTCATGTCTA\n+\nqaeiidaijsidaj\n'
        fpath_ill = os.path.join(original_reads_dir, 'pt_illumina.sfastq')
        reads_ill = '@seq2\nGCTACATGTCTA\n+\nqaeiidasidaj\n'
        open(fpath_454, 'w').write(reads_454)
        open(fpath_ill, 'w').write(reads_ill)

        #use only the 454 and sanger reads
        analsysis_config = {
            'inputs_filter':lambda input: input['platform'] in ('454',
                                                                'illumina')}
        do_analysis(project_settings=settings_path,
                    kind='clean_reads',
                    analysis_config=analsysis_config)

        do_analysis(project_settings=settings_path,
                    kind='prepare_mira_assembly',
                    analysis_config=analsysis_config)

        do_analysis(project_settings=settings_path,
                    kind='mira_assembly',
                    analysis_config=analsysis_config)


        test_dir.close()

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
