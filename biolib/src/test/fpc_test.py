'''
Created on 23/09/2009

@author: jose
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

import unittest, os

import biolib
from biolib.fpc import FPCMap

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class TestFPC(unittest.TestCase):
    'It tests the fpc functionality'

    @staticmethod
    def test_fpc():
        'It tests the fpc parsing'
        fpc_fname = os.path.join(DATA_DIR, 'biofpc.fpc')
        fpc = FPCMap(open(fpc_fname))
        assert fpc.name == 'demo'
        assert fpc.version == '8.5.1'

ok $mapio = Bio::MapIO->new(-format => 'fpc',
                                    -file   => Bio::Root::IO->catfile
                                    ('t','data','ctgdemo.fpc'));

$map = $mapio->next_map;

ok($map->length, 0);
ok($map->version, 7.2);
ok($map->modification_user, 'cari');
ok($map->group_type, 'Chromosome');
ok($map->group_abbr, 'Chr');
ok($map->core_exists, 0);

$count = 0;
foreach my $marker ($map->each_markerid) {
    $count++;
}

ok($count,150);

# add tests for get_markerobj

$count = 0;
foreach my $clone ($map->each_cloneid) {
    $count++;
}

ok($count,618);

# add tests for get_cloneobj

$count = 0;
foreach my $contig ($map->each_contigid) {
    $count++;
}

ok($count,2);

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestFPC.test_fpc']
    unittest.main()