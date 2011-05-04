'''
Created on 2009 eka 1

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

from franklin.db.databases_info import DatabasesInfo
import unittest
import os
from franklin.utils.misc_utils import TEST_DATA_DIR

class DatabasesInfoTest(unittest.TestCase):
    '''It tests the databases info utiliti '''
    @staticmethod
    def test_databases_info():
        '''It tests common use '''
        iabbrs = DatabasesInfo(os.path.join(TEST_DATA_DIR, 'GO.xrf_abbs'))
        assert iabbrs.dictionary.keys()
    @staticmethod
    def test_database_url_syntax():
        ''' It tests database url syntax'''
        db_info = DatabasesInfo(os.path.join(TEST_DATA_DIR, 'GO.xrf_abbs'))
        url = db_info.url_syntax('AGI_LocusCode', '123')
        good_url = "getEntry.pl?db_pick=[ChickGO/MaizeGO]&uid=123"
        assert good_url == url[38:]
        url = db_info.url_syntax('AGI_Loc', '123')
        assert url is None
    @staticmethod
    def test_last_db():
        '''It teste the last db in the file '''
        db_info = DatabasesInfo(os.path.join(TEST_DATA_DIR, 'GO.xrf_abbs'))
        assert db_info.dictionary['ASPGD_REF']

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
