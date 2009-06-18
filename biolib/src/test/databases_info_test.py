'''
Created on 2009 eka 1

@author: peio
'''
from biolib.db.databases_info import DatabasesInfo
import unittest
import os
import biolib
DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class DatabasesInfoTest(unittest.TestCase):
    '''It tests the databases info utiliti '''
    @staticmethod
    def test_databases_info():
        '''It tests common use '''
        iabbrs = DatabasesInfo(os.path.join(DATA_DIR, 'GO.xrf_abbs'))
        assert iabbrs.dictionary.keys()
    @staticmethod
    def test_database_url_syntax():
        ''' It tests database url syntax'''
        db_info = DatabasesInfo(os.path.join(DATA_DIR, 'GO.xrf_abbs'))
        url = db_info.url_syntax('AGI_LocusCode', '123')
        good_url = "getEntry.pl?db_pick=[ChickGO/MaizeGO]&uid=123"
        assert good_url == url[38:]
        url = db_info.url_syntax('AGI_Loc', '123')
        assert url is None
    @staticmethod
    def test_last_db():  
        '''It teste the last db in the file '''
        db_info = DatabasesInfo(os.path.join(DATA_DIR, 'GO.xrf_abbs'))
        assert db_info.dictionary['ASPGD_REF']
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()