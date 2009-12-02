'''
This module contain cleaners and filters to accept or reject founded seq
variations

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

import copy

from biolib.utils.seqio_utils import FileSequenceIndex
from biolib.snv.snv import (cap_enzymes, SNP, reference_variability,
                            get_reference_name, aggregate_alleles)
from biolib.seq.seq_analysis import infer_introns_for_cdna
from biolib.seq import seq_filters
from biolib.utils.cmd_utils import create_runner
from biolib.alignment_search_result import (FilteredAlignmentResults,
                                            get_alignment_parser)

SEQUENCE_FILTERS = {'aligner': seq_filters.create_aligner_filter,
                    'similar_seqs': seq_filters.create_similar_seqs_filter}

def create_reference_filter(seq_filter, filter_args):
    '''It filters the all snv that have a reference that is filtered by the
    given sequence filter'''
    seq_filter = SEQUENCE_FILTERS[seq_filter](**filter_args)
    cache = {}
    cache['seq_id'] = None
    cache['result'] = None
    def snv_filter(snv):
        'An snv filter that filters according to the reference sequence'
        if snv is None:
            return False
        snv = snv[0]
        #do we filter the reference sequence?
        sequence = snv.reference
        cache_id = sequence.id + sequence.name
        if cache_id == cache['seq_id']:
            return cache['result']
        else:
            result = seq_filter(sequence)
            cache['result'] = result
            cache['seq_id'] = cache_id
            return result
    return snv_filter

def create_unique_contiguous_region_filter(distance, genomic_db):
    '''It returns a filter that removes snv in a region that give more than one
    match or more than one match_parts'''
    parameters = {'database': genomic_db, 'program':'blastn'}
    blast_runner = create_runner(kind='blast', parameters=parameters)
    blast_parser = get_alignment_parser('blast')
    filters = [{'kind'           : 'min_scores',
            'score_key'      : 'similarity',
            'min_score_value': 90,
           },
           {'kind'           : 'min_length',
            'min_length_bp'  : 20,
           }
          ]
    def unique_contiguous_region_filter(snv):
        '''It filters out the snv in regions repeated in the genome or
        discontiguous'''
        snv = snv[0]
        #we make a blast
        #with the sequence around the snv
        location = snv.location
        start = location - distance
        end   = location + distance
        if start < 0:
            start = 0
        blast_fhand = blast_runner(snv.reference[start:end])[0]
        #now we parse the blast
        blast_result = blast_parser(blast_fhand)
        alignments = FilteredAlignmentResults(match_filters=filters,
                                          results=blast_result)
        #are there any similar sequences?
        try:
            alignment = alignments.next()
        except StopIteration:
            #if there is no similar sequence we assume that is unique
            return True
        #how many matches, it should be only one
        num_hits = len(alignment['matches'])
        if num_hits > 1:
            return False
        #how many match parts have the first match
        num_hsps = len(alignment['matches'][0]['match_parts'])
        if num_hsps == 1:
            return True
        else:
            return False
    return unique_contiguous_region_filter

def create_close_to_intron_filter(distance, genomic_db,
                                  genomic_seqs_fhand=None):
    '''It returns a filter that filters snv by the proximity to introns.

    The introns are predicted by blasting against a similar species'''
    introns_cache = {}
    introns_cache['seq_id'] = None
    introns_cache['introns'] = None

    genomic_seqs_index = None
    if genomic_seqs_fhand:
        genomic_seqs_index = FileSequenceIndex(genomic_seqs_fhand, 'fasta')

    def close_to_intron_filter(snv):
        'The filter'
        if snv is None:
            return False
        snv = snv[0]
        location = snv.location
        #where are the introns?
        sequence = snv.reference
        cache_id = sequence.id + sequence.name
        if cache_id == introns_cache['seq_id']:
            introns = introns_cache['introns']
        else:
            introns = infer_introns_for_cdna(sequence=snv.reference,
                                          genomic_db=genomic_db,
                                          genomic_seqs_index=genomic_seqs_index)
            introns_cache['introns'] = introns
            introns_cache['seq_id'] = cache_id
        for intron in introns:
            if abs(location - intron) < distance:
                return False
        return True
    return close_to_intron_filter

#filters
def create_high_variable_region_filter(max_variability, window=None):
    'It creates a filter that filters seq_vars from high variable regions'
    def high_variable_region_filter(snv):
        'The filter'
        if snv is None:
            return False
        (snv, context) = snv
        region_variability = reference_variability(snv, context, window=window)
        if region_variability > max_variability:
            return False
        else:
            return True
    return high_variable_region_filter

def create_close_to_seqvar_filter(distance):
    '''It returns a filter that filters snv by the proximity to other snvs.

    If the snv has another snv closer than DISTANCE, then this snv is
    filtered out'''
    def close_to_snv_filter(snv):
        'The filter'
        if snv is None:
            return False
        (snv, context) = snv
        location = snv.location
        for snv_in_contig in context:
            dist = abs(snv_in_contig.location - location)
            if dist == 0 :
                continue
            if dist < distance:
                return False
        return True
    return close_to_snv_filter

def create_snv_close_to_limit_filter(max_distance):
    '''This function is a function factory. These functions are filters than
     return true if those seqvars ar not clser to the limita than max_distance
     '''
    def snv_close_to_reference_limit(snv):
        '''True if the sequence variation is close to the contig limits.

        It checks if the seq_variation is not within the range covered by the
        consensus and if it's close to one of its limits. In both cases it will
        return True.
        '''
        if snv is None:
            return False
        (snv, context) = snv
        #where does the sequence variation starts and ends?
        location = snv.location

        #where does the consensus starts and ends?
        reference = snv.reference
        if reference  is None:
            msg = 'SeqVariation is not associated with an reference sequence'
            raise ValueError(msg)

        reference_len = len(reference)
        if ((reference_len - location  > max_distance) and
            (location  > max_distance)):
            return True
        else:
            return False

    return snv_close_to_reference_limit


#cleaners
def create_major_allele_freq_filter(frequency, libraries=None):
    '''It creates a cleaner in which the alleles from each library are removed
    if the most abundant allele frequency is bigger than the given one'''
    def major_allele_freq_cleaner(snv):
        'The cleaner'
        if snv is None:
            return False
        snv = snv[0]
        for library_info in snv.per_lib_info:
            if libraries and library_info['library'] not in libraries:
                continue
            alleles = library_info['alleles']
            read_number = 0
            first_read  = alleles[0]['reads']
            for allele in alleles:
                read_number += allele['reads']
            first_percent = (first_read) / float(read_number)
            if first_percent < frequency:
                return True
        return False

    return major_allele_freq_cleaner

def create_bad_quality_reads_cleaner(qual_treshold):
    '''This function is a factory function that giving a seqvar removes bad
    quality reads from it'''
    def bad_quality_reads_cleaner(snv):
        "The cleaner"
        if snv is None:
            return None
        (snv, context) = snv
        new_library_alleles = []
        for library_info in snv.per_lib_info:
            alleles = library_info['alleles']
            new_alleles = []
            for allele in alleles:
                new_qual = []
                if 'qualities' in allele:
                    for qual in allele['qualities']:
                        if qual > qual_treshold or qual is None:
                            new_qual.append(qual)
                else:
                    new_alleles.append(allele)
                    continue

                len_new_qual = len(new_qual)
                if len_new_qual == 0:
                    continue
                elif len_new_qual != len(allele['qualities']):
                    new_allele = {}
                    new_allele['allele']  = allele['allele']
                    new_allele['kind']    = allele['kind']
                    new_allele['reads']   = len(new_qual)
                    new_allele['qualities'] = new_qual
                    new_alleles.append(new_allele)
                else:
                    new_alleles.append(allele)

            # if we have not remove all the alleles, we add the allele to this
            # library
            if len(new_alleles) != 0:
                new_library = {}
                if 'library' in library_info:
                    new_library['library'] = library_info['library']
                new_library['alleles'] = new_alleles
                new_library_alleles.append(new_library)

        if new_library_alleles:
            return (snv.copy(per_lib_info=new_library_alleles), context)

    return bad_quality_reads_cleaner

def create_is_variable_in_some_filter(ignore_svns_respect_ref=True,
                                      libraries=None):
    '''This function factory creates a filter. This filter pass the svn if one
    of the given libraries is variable'''
    num_alleles = 2
    def is_variable_filter(snv):
        'The filter'
        if snv is None:
            return False
        snv = snv[0]
        for library_info in snv.per_lib_info:
            if libraries and library_info['library'] not in libraries:
                continue
            alleles = library_info['alleles']
            if ((not ignore_svns_respect_ref and len(alleles) == 1 and
                 alleles[0]['kind'] == SNP) or (len(alleles) >= num_alleles)):
                return True
        return False
    return is_variable_filter

def create_is_variable_in_aggregate_filter(libraries=None):
    '''This function factory creates a filter. This filter pass the svn if the
    aggregate produced by the alleles of all libraries is variable'''
    num_alleles = 2
    def is_variable_filter(snv):
        'The filter'
        if snv is None:
            return False
        snv = snv[0]
        alleles = []
        for library_info in snv.per_lib_info:
            if libraries and library_info['library'] not in libraries:
                continue
            alleles.extend(library_info['alleles'])
        alleles = aggregate_alleles(alleles)
        if len(alleles) >= num_alleles:
            return True
        return False
    return is_variable_filter

def create_kind_filter(kinds):
    'This filter filters by kind. It only passes snps or complex snv kinds'
    def kind_filter(snv):
        'The filter'
        if snv is None:
            return False
        snv = snv[0]
        if snv.kind in kinds:
            return True
        return False
    return kind_filter

def create_min_qual_per_lib_allele_cleaner(min_quality=70, default_quality=25):
    '''This function factory remove alleles that have not enough reads with good
    quality reads'''
    def qual_per_library_cleaner(snv):
        'The cleaner'
        if snv is None:
            return None
        (snv, context) = snv
        new_library_alleles = []
        for library_info in snv.per_lib_info:
            alleles = library_info['alleles']
            new_alleles = []
            for allele in alleles:
                qual_sum = 0
                for qual in allele['qualities']:
                    if qual is None:
                        qual = default_quality
                    qual_sum += qual
                if qual_sum  >= min_quality:
                    new_alleles.append(allele)

            # if we have not remove all the alleles, we add the allele to this
            # library
            if len(new_alleles) != 0:
                new_library = {}
                if 'library' in library_info:
                    new_library['library'] = library_info['library']
                new_library['alleles'] = new_alleles
                new_library_alleles.append(new_library)
        if new_library_alleles:
            return (snv.copy(per_lib_info=new_library_alleles), context)
        else:
            return None
    return qual_per_library_cleaner

def create_alleles_n_cleaner():
    '''This function factory removes alleles that composed by N'''
    def alleles_n_cleaner(snv):
        'The cleaner'
        if snv is None:
            return None
        (snv, context) = snv
        new_library_alleles = []
        for library_info in snv.per_lib_info:
            alleles = library_info['alleles']
            new_alleles = []
            for allele in alleles:
                if allele['allele'].lower() != 'n'*len(allele['allele']):
                    new_alleles.append(allele)

            # if we have not remove all the alleles, we add the allele to this
            # library
            if len(new_alleles) != 0:
                new_library = {}
                if 'library' in library_info:
                    new_library['library'] = library_info['library']
                new_library['alleles'] = new_alleles
                new_library_alleles.append(new_library)
        if new_library_alleles:
            return (snv.copy(per_lib_info=new_library_alleles), context)
    return alleles_n_cleaner

def create_cap_enzyme_filter(all_enzymes):
    '''This funtion is a factory function that creates a function that look
    if the seqvar is differently afected by some enzymes'''
    def enzymes_filter(snv):
        'The real filter'
        if snv is None:
            return False
        snv = snv[0]
        enzymes = cap_enzymes(snv, all_enzymes)
        if len(enzymes) != 0:
            return True
        else:
            return False
    return enzymes_filter

def create_reference_list_filter(references):
    'It creates a filter that removes the snvs with a reference not in the list'
    def reference_list_filter(snv):
        'It filters out the snvs with a reference not in the list'
        if snv is None:
            return False
        snv = snv[0]
        ref_name = get_reference_name(snv.reference)
        if ref_name in references:
            return True
        else:
            return False
    return reference_list_filter

def _allele_with_required_quality(allele, min_quality, default_quality):
    'It returns True if the allele has the mininum quality'
    #the quals with the default qualities
    quals = []
    for qual in allele['qualities']:
        if qual is None:
            qual = default_quality
        quals.append(qual)
    #we sort with the best qualities first
    qual_oris = zip(quals, allele['orientations'])
    qual_oris = sorted(qual_oris, lambda qo1, qo2: qo2[0] - qo1[0])
    #qualities for both orientations
    forward_quals, rev_quals = [], []
    for qual, ori in qual_oris:
        if ori:
            forward_quals.append(qual)
        else:
            rev_quals.append(qual)
    #now we calculate the qualities:
    qual = 0
    if forward_quals:
        qual += forward_quals[0]
        if len(forward_quals) > 1:
            qual += forward_quals[1] / 4.0
            if len(forward_quals) > 2:
                qual += forward_quals[2] / 4.0
    if rev_quals:
        qual += rev_quals[0]
        if len(rev_quals) > 1:
            qual += rev_quals[1] / 4.0
            if len(rev_quals) > 2:
                qual += rev_quals[2] / 4.0
    #the final check
    if qual >= min_quality:
        return True
    else:
        return False

def _alleles_with_required_quality(alleles, min_quality, default_quality):
    'It returns a list with the alleles have the minimum quality'
    ok_alleles = []
    for allele in alleles:
        if _allele_with_required_quality(allele, min_quality, default_quality):
            ok_alleles.append((allele['allele'], allele['kind']))
    return ok_alleles

def create_aggregate_allele_qual_cleaner(min_quality=45, default_quality=25):
    'It removes alleles with at least min quality.'
    def allele_qual_cleaner(snv):
        'It removes the alleles which have not enough quality'
        if snv is None:
            return None
        (snv, context) = snv
        new_library_alleles = []
        ag_alleles = snv.aggregate_alleles()
        #which alleles have the correct qualities
        ok_alleles = _alleles_with_required_quality(ag_alleles, min_quality,
                                                    default_quality)
        for library_info in snv.per_lib_info:
            alleles = library_info['alleles']
            new_alleles = []
            for allele in alleles:
                if (allele['allele'], allele['kind']) in ok_alleles:
                    new_alleles.append(allele)
            if new_alleles:
                new_library = copy.deepcopy(library_info)
                new_library['alleles'] = new_alleles
                new_library_alleles.append(new_library)
        if new_library_alleles:
            return (snv.copy(per_lib_info=new_library_alleles), context)
    return allele_qual_cleaner

