'''
Created on 15/01/2010

@author: peio
'''
# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

import os, unittest
from biolib.utils.misc_utils import DATA_DIR
from biolib.seq.seq_annotation import get_orthologs, get_descriptions_from_blasts

class OrthologsTests(unittest.TestCase):
    'It test basic use of the orthologs functions'
    @staticmethod
    def test_get_orthologs():
        'It tests the get_orthologs functions'
        blast_file  = open(os.path.join(DATA_DIR, 'melon_tair.xml'))
        blast_file2 = open(os.path.join(DATA_DIR, 'tair_melon.xml'))

        orthologs = get_orthologs(blast_file, blast_file2)
        #print orthologs.next()
        assert orthologs.next() == ('melon1', 'tair1')
        assert orthologs.next() == ('melon2', 'tair2')


class AnnotationTests(unittest.TestCase):
    'Annotations tests'
    @staticmethod
    def test_get_description_basic():
        'It tests if we can get description for seqs in blasts'
        # this fasta does not have definitio information
        blast = open(os.path.join(DATA_DIR, 'tair_melon.xml'))
        assert get_descriptions_from_blasts([blast]) == {}

    @staticmethod
    def test_get_description_with_funct():
        'It tests if we can get description for seqs in blasts. with mod funct'

        def mod_func(definition):
            'Arabidopsis def modifier'
            return definition.split('|')[2]
        # test with a modifier function
        blast_fhand    = open(os.path.join(DATA_DIR, 'blast2.xml'))
        blast = {'fhand':blast_fhand, 'desc_modifier': mod_func}
        assert get_descriptions_from_blasts([blast]) == \
                    {u'CUTC021854': u' ankyrin repeat family protein ',
                     u'CUTC021853': u' DNA-binding protein-related '}



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
