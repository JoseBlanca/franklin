'''
Steps to use with seqs in pipelines

Created on 03/12/2009

@author: peio
'''
from franklin.snv.snv_cleaner import (create_cap_enzyme_filter,
                                    create_snv_close_to_limit_filter,
                                    create_high_variable_region_filter,
                                    create_close_to_seqvar_filter,
                                    create_major_allele_freq_filter,
                                    create_is_variable_in_some_filter,
                                    create_bad_quality_reads_cleaner,
                                    create_min_qual_per_lib_allele_cleaner,
                                    create_alleles_n_cleaner,
                                    create_kind_filter,
                                    create_is_variable_in_aggregate_filter,
                                    create_reference_list_filter,
                                    create_aggregate_allele_qual_cleaner,
                                    create_close_to_intron_filter,
                                    create_reference_filter,
                                    create_unique_contiguous_region_filter)

#pylint:disable-msg=C0103
# Snp cleaning /filtering ####
snp_high_variable_region_filter = {
                         'function':create_high_variable_region_filter,
                         'arguments':{'max_variability':0.05},
                         'type':'filter',
                         'name':'high_variable_region',
                         'comment': 'It filters snp in a high variable regions'}

snp_close_to_seqvar_filter = {'function':create_close_to_seqvar_filter,
                            'arguments':{'distance':20},
                            'type':'filter',
                            'name':'close_to_seqvar',
                            'comment': 'It filters snp if it has a seqvar near'}


snp_major_allele_freq_filter = {'function':create_major_allele_freq_filter,
                                     'arguments':{'frequency':0.8},
                                     'type':'filter',
                                     'name':'major_allele_frec',
                              'comment':'It filters by mayor allele frequency'}

snp_cap_enzyme_filter  = {'function':  create_cap_enzyme_filter,
                      'arguments': {'all_enzymes':True},
                      'type':      'filter',
                      'name':      'enzyme_filter',
                      'comment':  'It filters by enzyme'}

snp_filter_is_variable_in_some = {'function': create_is_variable_in_some_filter,
                            'arguments':{},
                            'type':'filter',
                            'name':'variable_in_some',
                            'comment': 'It filters by allele_quantity'}

snp_filter_is_variable_in_aggregate = {
                          'function': create_is_variable_in_aggregate_filter,
                          'arguments':{},
                          'type':'filter',
                          'name':'variable_in_aggregate',
                          'comment': 'It filters by aggregate variability'}

snp_remove_by_read_number = {'function': create_min_qual_per_lib_allele_cleaner,
                            'arguments':{'num_reads':3},
                            'type':'mapper',
                            'name':'read_number',
                            'comment': 'It cleans alleles  by read number'}

snp_close_to_limit_filter = {'function':create_snv_close_to_limit_filter ,
                            'arguments':{'max_distance':12},
                            'type':'filter',
                            'name':'close_to_limit',
                            'comment': 'It filters by proximity to limit'}

snp_no_bad_quality_alleles_per_lib = \
                                 {'function':create_bad_quality_reads_cleaner,
                                  'arguments':{'qual_treshold':20},
                                  'type':'mapper',
                                  'name':'bad_quality_allele_striper',
                                  'comment': 'It removes bad quality reads'}
snp_no_baq_quality_alleles_agg = {'function':
                                          create_aggregate_allele_qual_cleaner,
                                  'arguments':{},
                                  'type':'mapper',
                                  'name':'bad_quality_allele_striper',
                                  'comment': 'It removes bad quality reads'}
snp_no_bad_quality_alleles_per_lib = {'function':
                                         create_min_qual_per_lib_allele_cleaner,
                                   'arguments':{},
                                  'type':'mapper',
                                  'name':'bad_quality_allele_striper',
                                  'comment': 'It removes bad quality reads'}
snp_remove_alleles_n = {'function':create_alleles_n_cleaner,
                        'arguments':{},
                        'type':'mapper',
                        'name':'allele_n_cleaner',
                        'comment': 'It removes alleles composed by N'
                        }
snp_filter_by_kind = {'function':create_kind_filter ,
                      'arguments':{},
                      'type':'filter',
                      'name':'kind_filter',
                      'comment': 'It filters by snv kind'}
snp_filter_reference_not_in_list = {
                        'function':create_reference_list_filter,
                        'arguments':{},
                        'type':'filter',
                        'name':'reference_list_filter',
                'comment': 'It filters by if the reference is not in the list'}
snp_filter_by_intron_proximity = {'function':create_close_to_intron_filter,
                               'arguments':{'distance':60, 'genomic_db':None,
                                            'genomic_seqs_fhand':None},
                               'type':'filter',
                               'name':'intron_filter',
                               'comment': 'It filters by proximity to a intron'}
msg  = 'It filters the snv in regions that give more than one hit or one hsp'
snv_filter_unique_contiguous_region = \
                            {'function':create_unique_contiguous_region_filter,
                           'arguments':{'distance':60, 'genomic_db':None,
                                        'genomic_seqs_fhand':None},
                           'type':'filter',
                           'name':'unique_contiguous_region',
                           'comment': msg}
snv_filter_by_reference = {'function':create_reference_filter,
                           'arguments':{'seq_filter':None, 'filter_args':None},
                           'type':'filter',
                           'name':'reference_filter',
                           'comment': 'It filters by snv.refernce properties'}

################################################################################
# PIPELINES
################################################################################

SNVPIPELINES = {'snp_basic'  : [snp_remove_alleles_n,
                                snp_no_baq_quality_alleles_agg,
                                snp_no_bad_quality_alleles_per_lib,
                                snp_filter_is_variable_in_some],
                }
