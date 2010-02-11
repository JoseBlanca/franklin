'''
Created on 07/10/2009

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


from franklin.statistics import create_distribution
from franklin.snv.snv import reference_variability, major_allele_frequency
from franklin.utils.collections_ import item_context_iter, FileCachedList


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
                                'xlabel':'reference variablity (snv / 100pb)',
                                'ylabel': 'Number of snvs'},
                        distrib_fhand=distrib_fhand,
                        plot_fhand=plot_fhand,
                        range_=range_, low_memory=True)

def calculate_snv_distrib(snvs, distrib_info, distrib_fhand=None,
                          plot_fhand=None, range_=None):
    '''It calculates the major allele frequency snv distribution.

    Number of snvs versus major allele frequency
    '''
    values = FileCachedList(float)
    for snv in snvs:
        value = distrib_info['function'](snv)
        if isinstance(value, list):
            values.extend(value)
        else:
            values.append(value)
    return create_distribution(values,
                               labels=distrib_info['labels'],
                               distrib_fhand=distrib_fhand,
                               plot_fhand=plot_fhand,
                               range_=range_, low_memory=True)

DISTRIBUTIONS = {
                 'ref_variability':{
                            'function':calculate_ref_variability_ditrib,
                            'snv_iter_kind':'snv_contexts'
                 },
                 'maf_distrib':{
                            'function':major_allele_frequency,
                            'snv_iter_kind':'snv',
                            'labels':{
                                'title':'Major allele frequency distribution',
                                'xlabel':'maf',
                                'ylabel': 'Number of snvs'
                            }
                },
}

def snv_distrib(snvs, kind, window=None, distrib_fhand=None, plot_fhand=None,
                range_=None):
    'It calculates one snv distribution of the given kind'
    distribution_orders = DISTRIBUTIONS[kind]
    if distribution_orders['snv_iter_kind'] == 'snv_contexts':
        snvs = item_context_iter(snvs)
        return distribution_orders['function'](snvs, window=window,
                                        distrib_fhand=distrib_fhand,
                                        plot_fhand=plot_fhand,
                                        range_=range_)
    else:
        return calculate_snv_distrib(snvs, distribution_orders,
                                        distrib_fhand=distrib_fhand,
                                        plot_fhand=plot_fhand,
                                        range_=range_)
