'''
Created on 27/10/2009

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

import unittest
from StringIO import StringIO
from tempfile import NamedTemporaryFile

from franklin.gmod.cmap import cmap_to_gff, cmap_to_mcf

class CmapTest(unittest.TestCase):
    'The cmap functionality'

    def test_cmap_gff(self):
        'cmap gff3 input creation'
        species = {'accession':'cmelo',
                   'full_name':'Cucumis melo',
                   'common_name':'melon',
                   'display_order':1}
        marker1 = {'type':'marker',
                   'name':'marker1',
                   'alias':'also_marker_1'}
        marker2 = {'type':'marker_',
                   'name':'marker2'}
        markers = {'marker1':marker1, 'marker2':marker2}
        feat_loc = [{'feature':'marker1', 'start':5, 'end':10},
                    {'feature':'marker2', 'start':15}]
        map1 = {'accession': 'gl1',
                'name': 'group1',
                'display_order':1,
                'feature_locations':feat_loc}
        map_set1 = {'species':'cmelo',
                   'accession':'ms1',
                   'name':'map set 1',
                   'short_name':'ms1',
                   'type':'Genetic',
                   'unit_modifier':0.01,
                   'maps' : [map1]}
        feat_loc2 = [{'feature':'marker1', 'start':1, 'end':5},
                    {'feature':'marker2', 'start':20}]
        map2 = {'accession': 'gl2',
                'name': 'group2',
                'display_order':2,
                'feature_locations':feat_loc2}
        map_set2 = {'species':'cmelo',
                   'accession':'ms2',
                   'name':'map set 2',
                   'short_name':'ms2',
                   'type':'Genetic',
                   'unit_modifier':0.01,
                   'maps' : [map2]}
        cmap = {'features': markers,
                'map_sets':[map_set1, map_set2],
                'species':{'cmelo': species}}
        fhand = NamedTemporaryFile()
        cmap_to_gff(cmap, fhand)
        correspondences = '''
##cmap_corr evidence_type_acc=ANB;ID1=marker1_0;map_set_acc1=ms1;map_acc1=gl1;ID2=marker1_1;map_set_acc2=ms2;map_acc2=gl2;
##cmap_corr evidence_type_acc=ANB;ID1=marker2_0;map_set_acc1=ms1;map_acc1=gl1;ID2=marker2_1;map_set_acc2=ms2;map_acc2=gl2;
'''
        result = open(fhand.name).read()
        assert correspondences in result
        assert "feature_accs=ms2_group2_marker1" in  result


        fhand = NamedTemporaryFile()
        cmap_to_mcf(cmap, fhand)
        result = open(fhand.name).read()
        assert "marker2\t20\tC02" in result
        assert "qtls\nmarker1\t7.5\t7.5\t7.5\t7.5\tC01" in result

        # Repeated markers allowed
        feat_loc2 = [{'feature':'marker1', 'start':1, 'end':5},
                    {'feature':'marker2', 'start':20},
                    {'feature':'marker2', 'start':20}]
        map2 = {'accession': 'gl2',
                'name': 'group2',
                'display_order':2,
                'feature_locations':feat_loc2}
        map_set2 = {'species':'cmelo',
                   'accession':'ms2',
                   'name':'map set 2',
                   'short_name':'ms2',
                   'type':'Genetic',
                   'unit_modifier':0.01,
                   'maps' : [map2]}
        cmap = {'features': markers,
                'map_sets':[map_set2],
                'species':{'cmelo': species}}

        cmap_to_gff(cmap, fhand)
        result = open(fhand.name).read()
        assert 'ID=marker2_1' in result

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
