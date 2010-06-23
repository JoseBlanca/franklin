'''
Created on 26/10/2009

@author: jose
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

import unittest, os
from StringIO import StringIO

from franklin.gff import (write_gff, features_in_gff, get_gff_header,
                          add_dbxref_to_feature)
from franklin.utils.misc_utils import DATA_DIR

class GffOutTest(unittest.TestCase):
    'It tests the gff writer'

    @staticmethod
    def test_simple_output():
        'We can write a simple gff3 file'
        feat1 = {'seqid': 'ctg123',
                 'type':  'gene',
                 'start': 1000,
                 'end':   9000,
                 #'id':    'gene00001',
                 #'name':  'EDEN',
                 'attributes':{'ID':'gene00001', 'Name':'EDEN'}
                 }
        feats = ['##sequence-region ctg123 1 1497228', feat1]
        result = '''##gff-version 3
##sequence-region ctg123 1 1497228
ctg123\t.\tgene\t1000\t9000\t.\t.\t.\tID=gene00001;Name=EDEN\n'''
        outh = StringIO()
        write_gff(feats, outh)
        assert result == outh.getvalue()

        feat1 = {'id':'23',
                 'seqid': 'ctg123',
                 'type':  'gene',
                 'start': 1000,
                 'end':   9000,
                 'attributes' : {'Parent': ['p1', 'p2']}}
        feats = [feat1]
        result = '''##gff-version 3
ctg123\t.\tgene\t1000\t9000\t.\t.\t.\tID=23;Parent=p1,p2\n'''
        outh = StringIO()
        write_gff(feats, outh)
        assert result == outh.getvalue()

    @staticmethod
    def test_features_in_gff():
        'It test the features_in_gff function'
        gff_fhand = open(os.path.join(DATA_DIR, 'map_fis.gff3'))
        features  = list(features_in_gff(gff_fhand, 3))
        assert len(features) == 99
        assert features[1]['name'] == 'ctg0'
        assert features[98]['name'] == 'Cm13_B04'

    @staticmethod
    def test_get_gff_header():
        'It test the features_in_gff function'
        gff_fhand = open(os.path.join(DATA_DIR, 'map_fis.gff3'))
        header = get_gff_header(gff_fhand)
        assert len(header) == 1
        assert header[0] == '##gff-version 3'

    @staticmethod
    def test_add_():
        'It test the features_in_gff function'
        feature = {'seqid': 'ctg123',
                   'type':  'gene',
                   'start': 1000,
                   'end':   9000,
                   'attributes':{'ID':'gene00001', 'Name':'EDEN'}
                   }

        dbxref_db = 'test'
        dbxref_id = 'id100'
        add_dbxref_to_feature(feature, dbxref_db, dbxref_id)
        assert feature['attributes']['Dbxref'] == 'test:id100'

        feature = {'seqid': 'ctg123',
                   'type':  'gene',
                   'start': 1000,
                   'end':   9000,
                   'attributes':{'ID':'gene00001', 'Name':'EDEN',
                                 'Dbxref':'test2:id101'}
                   }
        add_dbxref_to_feature(feature, dbxref_db, dbxref_id)
        assert feature['attributes']['Dbxref'] == 'test2:id101,test:id100'


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
