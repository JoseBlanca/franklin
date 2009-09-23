'''
Created on 22/09/2009

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

from biolib.collections_ import FileCachedList
from biolib.statistics import create_distribution

def calculate_read_coverage(pileup, distrib_fhand=None, plot_fhand=None,
                            range_=None):
    '''Given a sam pileup file it returns the coverage distribution.

    The coverage shows how many times the bases has been read.
    '''
    coverages = FileCachedList(int)
    for line in pileup:
        if line.isspace():
            continue
        position_cov = line.split()[3]
        coverages.append(position_cov)
    #now the distribution
    return create_distribution(coverages,
                               labels={'title':'Read coverage distribution',
                                      'xlabel':'coverage',
                                      'ylabel': 'Number of positions'},
                               distrib_fhand=distrib_fhand,
                               plot_fhand=plot_fhand,
                               range_=range_, low_memory=True)
