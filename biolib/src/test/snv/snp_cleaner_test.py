'''
Created on 2009 uzt 30

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

import unittest
from StringIO import StringIO

from biolib.seqvar.seqvariation import (SNP, DELETION, INVARIANT, Snv, COMPLEX)
from biolib.seqvar.snp_cleaner import (#create_major_allele_freq_filter,
                                       create_close_to_seqvar_filter,
                                       create_cap_enzyme_filter,
                                       create_high_variable_region_filter,
                                       create_is_variable_in_some_filter,
                                       create_is_variable_in_aggregate_filter,
                                       create_bad_quality_reads_cleaner,
                                       create_snv_close_to_limit_filter,
                                       create_major_allele_freq_filter,
                                       create_read_number_cleaner,
                                       create_alleles_n_cleaner,
                                       create_kind_filter,
                                       create_reference_list_filter,
                                       create_allele_qual_cleaner)
from biolib.seqvar.sam_pileup import snv_contexts_in_sam_pileup
from biolib.pipelines import (pipeline_runner, snp_filter_is_variable_in_some,
                              snp_filter_is_variable_in_aggregate,
                              snp_filter_by_kind,
                              snp_filter_reference_not_in_list)
from biolib.seqs import SeqWithQuality

class SeqVariationFilteringTest(unittest.TestCase):
    'It checks the filtering methods.'

    @staticmethod
    def test_enzyme_filter():
        'It test the enzyme filter'
        seq  = 'ATGATGATG' + 'gaaattc' + 'ATGATGATGTGGGAT'
        reference = SeqWithQuality(seq=seq, name='ref')
        snp = Snv(per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':4, 'kind':INVARIANT},
                                  {'allele':'A', 'reads':2, 'kind':DELETION}]}],
                        location=11, reference=reference)
        cap_enzime = create_cap_enzyme_filter(True)
        snp = (snp, 'context')
        assert cap_enzime(snp) == True

    @staticmethod
    def test_major_allele_freq_filter_snv():
        'It test the first allele percent filter'
        snp = Snv(per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':4, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':2, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}],
                  location=11, reference='reference')
        snp = (snp, 'context')
        frec_cleaner = create_major_allele_freq_filter(0.51)
        assert frec_cleaner(snp)

        snp = Snv(location=11, reference='reference', per_lib_info=[
                        {'library':'somelib',
                         'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'library':'somelib2',
                         'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])
        snp = (snp, 'context')
        frec_cleaner = create_major_allele_freq_filter(0.51)
        assert not frec_cleaner(snp)

    @staticmethod
    def test_high_variable_region_filter():
        'It test percent_variations_in_seq_ref filter'
        reference = SeqWithQuality(seq='atatat')
        snp = Snv(location=1, reference=reference, per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp1 = Snv(location=4, reference=reference, per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp2 = Snv(location=7, reference=reference, per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        context = [snp, snp1, snp2]
        snp     = (snp1, context)
        # This reference variability is 50
        filter40 = create_high_variable_region_filter(40)
        assert not filter40(snp)

        filter60 = create_high_variable_region_filter(60)
        assert filter60(snp)

    @staticmethod
    def test_create_filter_close_seqvar():
        'It test the first allele percent filter'
        reference = 'reference'
        snp = Snv(location=10, reference=reference, per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp1 = Snv(location=15, reference=reference, per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp2 = Snv(location=20, reference=reference, per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])
        context = [snp, snp1, snp2]
        snp = (snp, context)

        filter4 = create_close_to_seqvar_filter(5)
        filter7 = create_close_to_seqvar_filter(7)
        assert filter4(snp)
        assert not filter7(snp)

    @staticmethod
    def test_allele_number_filter():
        'It test percent_variations_in_seq_ref filter'
        reference = 'atatat'
        snp = Snv(location=10, reference=reference, per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])



        snp1 = Snv(location=15, reference=reference, per_lib_info=[
                    {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT}]}])

        snp2 = Snv(location=20, reference=reference, per_lib_info=[
                      {'library':'lib1',
                       'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                  {'allele':'T', 'reads':3, 'kind':SNP}]},
                      {'library':'lib2',
                       'alleles':[{'allele':'T', 'reads':2,'kind':SNP}]}])

        allele_number_filter = create_is_variable_in_some_filter()
        assert  allele_number_filter((snp, [snp]))
        assert not allele_number_filter((snp1, [snp]))
        assert allele_number_filter((snp2, [snp]))

        allele_number_filter = create_is_variable_in_some_filter(
                                                           libraries=['lib2'])
        assert not allele_number_filter((snp2, [snp]))

    @staticmethod
    def test_is_variable_in_aggregate_filter():
        'It test percent_variations_in_seq_ref filter'
        reference = 'atatat'
        snv = Snv(location=10, reference=reference, per_lib_info=[
                    {'library':'lib1',
                     'alleles':[{'allele':'A', 'reads':3, 'kind':SNP,
                                 'qualities':[]}]},
                    {'library':'lib2',
                     'alleles':[{'allele':'T', 'reads':2, 'kind':SNP,
                                 'qualities':[]}]},
                    {'library':'lib3',
                     'alleles':[{'allele':'T', 'reads':2, 'kind':SNP,
                                 'qualities':[]}]}])
        snv = snv, [snv]
        filter_ = create_is_variable_in_aggregate_filter()
        assert filter_(snv)

        filter_ = create_is_variable_in_aggregate_filter(libraries=['lib3',
                                                                    'lib2'])
        assert not filter_(snv)


    @staticmethod
    def test_create_seqvar_close_to_limt():
        'It tests the close to limit filter function '
        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])
        snp = (snp, 'context')
        filter_ = create_snv_close_to_limit_filter(3)
        assert filter_(snp)
        filter_ = create_snv_close_to_limit_filter(5)
        assert not filter_(snp)

    @staticmethod
    def test_remove_bad_quality_reads():
        'It test it we can remove the alleles with bad qualities'
        #the good reads in 454 are around a quality of 25 to 40
        #a bad letter is less than 20
        #for sanger the phred quality values are related to the probability
        #of having an error in the sequence
        #http://www.phrap.com/phred/
        #Phred quality score  Probability base is wrong   Accuracy of the base
        #       10                   1 in 10                  90%
        #       20                   1 in 100                 99%
        #       30                   1 in 1,000               99.9%
        #       40                   1 in 10,000              99.99%
        #       50                   1 in 100,000             99.999%
        bad_quality_cleaner = create_bad_quality_reads_cleaner(25)
        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':4, 'kind':INVARIANT,
                                     'qualities':[30, 30, 30, 20]},
                                    {'allele':'T', 'reads':2,'kind':SNP,
                                     'qualities':[30, 20]}]}])
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]
        assert len(snp.per_lib_info[0]['alleles']) == 2
        assert snp.per_lib_info[0]['alleles'][0]['reads'] == 3

        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT,
                             'qualities':[30, 30, 30, 20]},
                            {'allele':'T', 'reads':2,'kind':SNP,
                             'qualities':[20, 20]}]}])
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]

        assert len(snp.per_lib_info[0]['alleles']) == 1
        assert snp.per_lib_info[0]['alleles'][0]['reads'] == 3

        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                {'alleles':[{'allele':'AT', 'reads':4, 'kind':DELETION,
                             'qualities':[None, None, None, None]},
                            {'allele':'T', 'reads':2,'kind':SNP,
                             'qualities':[30, 30]}]}])
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]

        assert len(snp.per_lib_info[0]['alleles']) == 2
        assert snp.per_lib_info[0]['alleles'][0]['reads'] == 4
        assert len(snp.per_lib_info[0]['alleles'][0]['qualities']) == 4

        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                {'alleles':[{'allele':'AT', 'reads':4, 'kind':DELETION},
                            {'allele':'T', 'reads':2,'kind':SNP }]}])
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]
        assert len(snp.per_lib_info[0]['alleles']) == 2



    @staticmethod
    def test_read_number_cleaner():
        'Tests read number_cleaner'
        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':4, 'kind':INVARIANT,
                                     'quality':[30, 30, 30, 20]},
                                    {'allele':'T', 'reads':1,'kind':SNP,
                                     'quality':[30, 20]}]}])
        snp = (snp, 'context')
        read_number_cleaner = create_read_number_cleaner(2)
        snp = read_number_cleaner(snp)[0]
        assert len(snp.per_lib_info[0]['alleles']) == 1

        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':1, 'kind':INVARIANT,
                                     'quality':[30, 30, 30, 20]},
                                    {'allele':'T', 'reads':1,'kind':SNP,
                                     'quality':[30, 20]}]}])
        snp = (snp, 'context')
        read_number_cleaner = create_read_number_cleaner(2)
        snp = read_number_cleaner(snp)
        assert snp is None

    @staticmethod
    def test_alleles_n_cleaner():
        'It checks the removal of alleles N or NN'
        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'N', 'reads':10},
                                    {'allele':'nN', 'reads':10},
                                    {'allele':'T', 'reads':10}]}])
        snp = (snp, 'context')
        n_cleaner = create_alleles_n_cleaner()
        snp = n_cleaner(snp)[0]
        assert len(snp.per_lib_info[0]['alleles']) == 1

    @staticmethod
    def test_alleles_quality_cleaner():
        'It test the new allele quality cleaner'
        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':2, 'kind':SNP,
                                     'qualities':[29, 22],
                                     'orientations':[True, True]}]}])
        filter_ = create_allele_qual_cleaner()
        snp = filter_((snp, [snp]))[0]
        assert len(snp.per_lib_info[0]['alleles']) == 1

        filter_ = create_allele_qual_cleaner(min_quality=32)
        assert filter_((snp, [snp])) is None

        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':2, 'kind':SNP,
                                     'qualities':[20, 22],
                                     'orientations':[True, False]}]}])
        filter_ = create_allele_qual_cleaner()
        snp = filter_((snp, [snp]))[0]
        assert len(snp.per_lib_info[0]['alleles']) == 1

        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':SNP,
                                     'qualities':[20, 22, None],
                                     'orientations':[True, True, False]}]}])
        filter_ = create_allele_qual_cleaner()
        snp = filter_((snp, [snp]))[0]
        assert len(snp.per_lib_info[0]['alleles']) == 1

        #we take into account the alleles from differnt libraries
        snp = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':1, 'kind':SNP,
                                     'qualities':[25],
                                     'orientations':[True]},
                                     {'allele':'T', 'reads':1, 'kind':SNP,
                                     'qualities':[25],
                                     'orientations':[True]},]},
                        {'alleles':[{'allele':'A', 'reads':1, 'kind':SNP,
                                     'qualities':[25],
                                     'orientations':[False]},
                                     {'allele':'T', 'reads':1, 'kind':SNP,
                                     'qualities':[25],
                                     'orientations':[True]}]}])
        filter_ = create_allele_qual_cleaner()
        snp = filter_((snp, [snp]))[0]
        assert len(snp.per_lib_info[0]['alleles']) == 1

    @staticmethod
    def test_kind_filter():
        'It checks the filter by snv kind'
        snv = Snv(location=4, reference='atatatatat', per_lib_info=[
                        {'alleles':[{'allele':'A', 'reads':10, 'kind':SNP},
                                    {'allele':'T', 'reads':10,
                                     'kind':INVARIANT}]}])
        snv = (snv, 'context')
        kind_filter = create_kind_filter([SNP, COMPLEX])
        assert kind_filter(snv)

    @staticmethod
    def test_reference_in_ref_list():
        'It can filter out the snvs with a reference not in the given list'
        allowed_refs = ['ref1']
        snv1 = Snv(location=4, reference='ref1', per_lib_info=[])
        snv2 = Snv(location=4, reference='ref2', per_lib_info=[])
        filter_ = create_reference_list_filter(allowed_refs)

        assert filter_((snv1, [snv1]))
        assert not filter_((snv2, [snv2]))

    @staticmethod
    def xtest_svn_pipeline():
        'The complete snv minning process from pileup to filtered Snvs'
        tpileup1 = '''ref1     1      A      8       NTtGGaCC       ~==<<~>>
ref1     2      A      1       ,,,,,,,       aaaaaaa
ref1     3      A      1       ,       a
ref1     4      A      1       ,       a'''
        tpileup2 = '''ref1     1      A      3       TTt       ccc
ref1     2      A      1       ,,,,,,,       aaaaaaa
ref1     3      A      1       ,       a
ref1     4      A      1       ,       a'''
        tpileup3 = '''ref1     1      A      3       GGG       ddd
ref1     2      A      1       ,,,,,,,       aaaaaaa
ref1     3      A      1       ,       a
ref1     4      A      1       ,       a'''

        pileup1 = StringIO(tpileup1)
        pileup2 = StringIO(tpileup2)
        pileup3 = StringIO(tpileup3)

        #we get the snvs from the pileups
        snv_contexts = list(snv_contexts_in_sam_pileup([pileup3,
                                                        pileup2, pileup1],
                                        libraries=['lib3', 'lib2', 'lib1']))

        #we filter them
        pipeline = []
        config = {}
        # variable in lib1
        pipeline.append(snp_filter_is_variable_in_some)
        config['variable_in_some'] = {'libraries':'lib1'}

        # variable in the aggregate of all other libraries
        pipeline.append(snp_filter_is_variable_in_aggregate)
        config['variable_in_aggregate'] = {'libraries':['lib2', 'lib3']}
        # kind snp
        pipeline.append(snp_filter_by_kind)
        config['kind_filter'] = {'kinds':[SNP]}
        # variability of the unigene
        #pipeline.append(snp_high_variable_region_filter)
        #config['high_variable_region'] = {'max_variability':3}
        # close to seqvar filter
        #pipeline.append(snp_close_to_seqvar_filter)
        #config['close_to_seqvar'] = {'distance':30}
        #close to limit filter
        #pipeline.append(snp_close_to_limit_filter)
        #config['close_to_limit'] = {'max_distance': 30}
        #cap filter
        #pipeline.append(snp_cap_enzyme_filter)
        #config['enzyme_filter'] = {'all_enzymes': False}
        #reference filter
        pipeline.append(snp_filter_reference_not_in_list)
        config['reference_list_filter'] = {'references':['ref1']}

        snv_contexts = pipeline_runner(pipeline, snv_contexts, config)

        snv_contexts = list(snv_contexts)
        #in lib3 the G is read just twice
        assert snv_contexts[0][0].per_lib_info[0]['alleles'][0]['reads'] == 3
        assert (snv_contexts[0][0].per_lib_info[0]['alleles'][0]['qualities'] ==
                                                                   [67, 67, 67])

        #the same, but with the snp basic pipeline
        pileup1 = StringIO(tpileup1)
        pileup2 = StringIO(tpileup2)
        pileup3 = StringIO(tpileup3)
        #we get the snvs from the pileups
        snv_contexts = list(snv_contexts_in_sam_pileup([pileup3,
                                                        pileup2, pileup1],
                                        libraries=['lib3', 'lib2', 'lib1']))
        snv_contexts = pipeline_runner('snp_basic', snv_contexts)
        snv_contexts = list(snv_contexts)
        #in lib3 the G is read just twice
        snv0 = snv_contexts[0][0]
        assert len(snv0.per_lib_info) == 3
        assert len(snv0.per_lib_info[2]['alleles']) == 2
        assert snv0.per_lib_info[2]['alleles'][0]['allele'] == 'C'
        assert snv0.per_lib_info[2]['alleles'][1]['allele'] == 'T'
        assert snv0.per_lib_info[2]['alleles'][0]['qualities'] == [29, 29]


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SeqVariationFilteringTest.test_svn_pipeline']
    unittest.main()
