'''
Created on 07/10/2009

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


from biolib.statistics import create_distribution
from biolib.seqvar.seqvariation import (reference_variability,
                                        major_allele_frequency)
from biolib.collections_ import item_context_iter, FileCachedList

def calculate_ref_variability_ditrib(snv_contexts, window=None,
                                     distrib_fhand=None, plot_fhand=None,
                                     range_=None):
    '''It calculates an snv distribution.

    How many snv are located in sequence reference with some variability.
    The variability is calculated as number of snps / lenght(context)
    '''
    variabilities = FileCachedList(float)
    for snv, context in snv_contexts:
        region_variability = reference_variability(snv, context, window=window)
        variabilities.append(region_variability)
    return create_distribution(variabilities,
                       labels={'title':'Snv reference variability distribution',
                                'xlabel':'variablity',
                                'ylabel': 'Number of snvs'},
                        distrib_fhand=distrib_fhand,
                        plot_fhand=plot_fhand,
                        range_=range_, low_memory=True)

def calculate_maf_distrib(snvs, distrib_fhand=None, plot_fhand=None,
                          range_=None):
    '''It calculates the major allele frequency snv distribution.

    Number of snvs versus major allele frequency
    '''
    mafs = FileCachedList(float)
    for snv in snvs:
        maf = major_allele_frequency(snv)
        mafs.extend(maf)
    return create_distribution(mafs,
                       labels={'title':'Major allele frequency distribution',
                                'xlabel':'maf',
                                'ylabel': 'Number of snvs'},
                        distrib_fhand=distrib_fhand,
                        plot_fhand=plot_fhand,
                        range_=range_, low_memory=True)
