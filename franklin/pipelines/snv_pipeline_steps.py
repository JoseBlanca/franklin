'''
Steps to use with seqs in pipelines

Created on 03/12/2009

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

from franklin.snv.snv_annotation import create_snv_annotator
from franklin.snv.snv_filters import (create_high_variable_region_filter,
                                      create_unique_contiguous_region_filter,
                                      create_close_to_intron_filter,
                                      create_close_to_snv_filter,
                                      create_snv_close_to_limit_filter,
                                      create_major_allele_freq_filter,
                                      create_kind_filter,
                                      create_cap_enzyme_filter,
                                      create_is_variable_filter,
                                      create_reference_in_list_filter,
                                      create_not_variable_in_group_filter,
                                      create_min_groups_filter,
    create_in_segment_filter)

snv_bam_annotator = {'function':create_snv_annotator,
          'arguments':{'bam_fhand':None, 'min_quality':45,
                       'default_sanger_quality':25,
                       'min_mapq':15,
                       'min_num_alleles':1},
          'type':'mapper',
          'name':'snv_bam_annotator',
          'comment': 'It annotates the snvs from a bam file'}

unique_contiguous_region_filter = {
          'function':create_unique_contiguous_region_filter,
          'arguments':{'distance':60, 'genomic_db':None,
                        'genomic_seqs_fpath':None},
          'type':'filter',
          'name':'uniq_contiguous',
          'comment': 'A blast in the near region gave several matches'}

close_to_intron_filter = {'function':create_close_to_intron_filter,
          'arguments':{'distance':60},
          'type':'filter',
          'name':'close_to_intron',
          'comment': 'An intron is located closer than N base pairs'}

high_variable_region_filter = {
          'function':create_high_variable_region_filter,
          'arguments':{'max_variability':0.05},
          'type':'filter',
          'name':'high_variable_region',
          'comment': 'The snv is in a region with more than N % of variability'}

close_to_snv_filter = {
          'function':create_close_to_snv_filter,
          'arguments':{'distance':60},
          'type':'filter',
          'name':'close_to_snv',
          'comment': 'The snv is closer than N nucleotides to another snv'}

close_to_limit_filter = {
          'function':create_snv_close_to_limit_filter ,
          'arguments':{'distance':12},
          'type':'filter',
          'name':'close_to_limit',
          'comment':'The snv is closer than N nucleotides to sequence limit'}

major_allele_freq_filter = {
          'function':create_major_allele_freq_filter,
          'arguments':{'frequency':0.8, 'groups':None, 'group_kind':None},
          'type':'filter',
          'name':'maf',
          'comment':'The more frequent allele is more frequent than frec'}

kind_filter = {
          'function':create_kind_filter,
          'arguments':{},
          'type':'filter',
          'name':'by_kind',
          'comment': 'It filters by snv kind'}

cap_enzyme_filter = {
          'function':  create_cap_enzyme_filter,
          'arguments': {'all_enzymes':True},
          'type': 'filter',
          'name': 'cap_enzyme',
        'comment': 'It filters by enzymes that recognize different snp alleles'}

is_variable_filter = {
          'function': create_is_variable_filter,
          'arguments': {'group_kind': None, 'groups':None},
          'type': 'filter',
          'name': 'is_variable',
          'comment': 'It filters by variability is selected groups'}

is_not_variable_filter = {
          'function': create_not_variable_in_group_filter,
          'arguments': {'group_kind': None, 'groups':None},
          'type': 'filter',
          'name': 'is_not_variable',
          'comment': 'It filters by not variables in selected groups'}
ref_not_in_list = {
          'function': create_reference_in_list_filter,
          'arguments': {'seq_list': None},
          'type': 'filter',
          'name': 'ref_not_in_list',
          'comment': 'Filters by given list of seq names'}
min_groups = {
          'function': create_min_groups_filter,
          'arguments': {'min_groups': None},
          'type': 'filter',
          'name': 'min_groups',
          'comment': 'Filters by the number of groups that support the SNV'}
in_segment = {
              'function': create_in_segment_filter,
              'arguments': {'segments': None},
              'type': 'filter',
              'name': 'in_segment',
              'comment': 'Filters snvs if they are inside the given segments '
              }
################################################################################
# PIPELINES
################################################################################

SNV_PIPELINES = {'snv_bam_annotator': [snv_bam_annotator]}

SNV_STEPS = [snv_bam_annotator, unique_contiguous_region_filter,
             close_to_intron_filter, high_variable_region_filter,
             close_to_snv_filter, close_to_limit_filter,
             major_allele_freq_filter, kind_filter, cap_enzyme_filter,
             is_variable_filter, is_not_variable_filter, ref_not_in_list,
             min_groups, in_segment]
