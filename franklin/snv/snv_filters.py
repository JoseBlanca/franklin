'''
Created on 19/02/2010

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

from Bio import SeqIO

from franklin.utils.cmd_utils import create_runner
from franklin.alignment_search_result import (FilteredAlignmentResults,
                                              get_alignment_parser)
from franklin.seq.seq_analysis import (infer_introns_for_cdna,
                                       similar_sequences_for_blast)
from franklin.seq.readers import guess_seq_file_format
from franklin.snv.snv_annotation import (calculate_maf_frequency,
                                         snvs_in_window, calculate_snv_kind,
                                         calculate_cap_enzymes,
                                         variable_in_groupping,
                                         _get_group)
from franklin.seq.seqs import get_seq_name

# In a filter TRUE result means that a snv does NOT pass the filter.
# So it writes it to the vcf

FILTER_DESCRIPTIONS = {
    'uniq_contiguous':
        {'id': 'UCR',
         'description':'Region is not unique or non contiguous'},
    'close_to_intron':
        {'id': 'I%2d',
         'description':'An intron is located closer than %2d base pairs'},
    'high_variable_reg':
        {'id': 'HVR%d',
    'description':'The region has more than %d snvs per 100 bases'},
    'close_to_snv':
        {'id':'cs%2d',
         'description':'The snv is closer than %d nucleotides to another snv'},
    'close_to_limit':
        {'id':'cl%2d',
      'description':'The snv is closer than %d nucleotides the reference edge'},
    'maf':
        {'id':'maf%.2f',
       'description':'The most frequent allele frequency is greater than %.2f'},
    'by_kind':
        {'id':'vk%s',
         'description':'It is not an %s'},
    'cap_enzymes':
        {'id':'ce%s',
         'description':'SNV is not a CAP detectable by the enzymes: %s'},
    'is_variable':
        {'id':'v%s%i',
        'description':'It is not variable in the %s : %s. All together: %s'},
    'is_not_variable':
        {'id':'nv%s%i',
        'description':'It is variable in the %s : %s. All together: %s'},
    'ref_not_in_list':
        {'id':'rnl',
        'description':'Filters by given list of seq names'},
    'min_groups':
        {'id':'m%s%i',
        'description':'SNV read in less than %i %s'},
    }

FILTER_COUNTS = {}

def get_filter_description(filter_name, parameters, filter_descriptions):
    'It returns the short id and the description'
    if (filter_name, parameters) in filter_descriptions:
        return filter_descriptions[filter_name, parameters]
    id_ = FILTER_DESCRIPTIONS[filter_name]['id']
    desc = FILTER_DESCRIPTIONS[filter_name]['description']

    if filter_name == 'by_kind':
        short_name, description = _get_nd_kind(id_, desc, parameters)
    elif filter_name == 'cap_enzymes':
        short_name, description = _get_nd_ce(id_, desc, parameters)
    elif filter_name == 'is_variable':
        short_name, description = _get_nd_iv(id_, desc, parameters)
    elif filter_name == 'is_not_variable':
        short_name, description = _get_nd_iv(id_, desc, parameters)
    elif filter_name == 'high_variable_reg':
        short_name, description = _get_nd_hvr(id_, desc, parameters)
    elif filter_name == 'min_groups':
        short_name, description = _get_min_groups_desc(id_, desc, parameters)
    else:
        if '%' in id_:
            short_name = id_ % parameters
        else:
            short_name = id_
        if '%' in desc:
            description = desc % parameters
        else:
            description = desc

    filter_descriptions[filter_name, parameters] = short_name, description

    return short_name, description

def _get_min_groups_desc(id_, desc, parameters):
    'It returns the name and id of the snv filter for min_groups'
    group_letter = parameters[1][0]
    min_group_num = parameters[0]
    short_name = id_ % (group_letter, min_group_num)
    description = desc % (min_group_num, parameters[1])

    return short_name, description

def _get_nd_hvr(id_, desc, parameters):
    'It returns the name and id of the snv filter for by is_variable filter'
    number = int(parameters[0] * 100)
    short_name = id_ % number
    description = desc % number

    return short_name, description

def _get_nd_iv(id_, desc, parameters):
    'It returns the name and id of the snv filter for by is_variable filter'
    global FILTER_COUNTS
    if desc not in FILTER_COUNTS:
        FILTER_COUNTS[desc] = 0
    FILTER_COUNTS[desc] += 1
    groups = {'libraries':'lb', 'read_groups':'rg', 'samples':'sm'}
    short_name = id_ % (groups[parameters[0]], FILTER_COUNTS[desc])
    groups = ','.join(parameters[1])
    description = desc % (parameters[0], groups, parameters[2])

    return short_name, description

def _get_nd_kind(id_, desc, parameters):
    'It returns the name and id of the snv filter for by kind filter'
    vkinds = {0:'snp', 1:'insertion', 2:'deletion', 3:'invariant', 4:'indel',
              5:'complex'}
    kind = vkinds[parameters]
    short_name = id_ % kind[0]
    description = desc % kind
    return short_name, description

def _get_nd_ce(id_, desc, parameters):
    'It returns the name and id of the snv filter for cap_enzyme filter'
    if parameters:
        enzymes = 'all'
        booltag = 't'
    else:
        enzymes = 'cheap ones'
        booltag = 'f'

    short_name = id_ % booltag
    description = desc % enzymes

    return short_name, description

def _add_filter_result(snv, filter_name, result, threshold=None):
    'It adds the filter to the SeqFeature qualifiers'
    qualifiers = snv.qualifiers
    if 'filters' not in qualifiers:
        qualifiers['filters'] = {}
    if filter_name not in qualifiers['filters']:
        qualifiers['filters'][filter_name] = {}
    qualifiers['filters'][filter_name][threshold] = result

def _get_filter_result(snv, filter_name, threshold=None):
    'It gets the result of a filter. Returns None if the filter is not done'
    qualifiers = snv.qualifiers
    if 'filters' not in qualifiers:
        return None
    try:
        result = qualifiers['filters'][filter_name][threshold]
        return result
    except KeyError:
        return None

def create_reference_in_list_filter(seq_list):
    '''It filters sequences looking in a list. If the sequence is in a list it
    passes the filter'''
    def reference_in_list_filter(sequence):
        "The filter"
        if sequence is None:
            return None

        name = get_seq_name(sequence)
        if name in seq_list:
            result = True
        else:
            result = False

        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'ref_not_in_list')
            if previous_result is not None:
                continue
            _add_filter_result(snv, 'ref_not_in_list', result)
        return sequence
    return reference_in_list_filter


def create_unique_contiguous_region_filter(distance, genomic_db,
                                           genomic_seqs_fpath):
    '''It returns a filter that removes snv in a region that give more than one
    match or more than one match_parts'''
    parameters = {'database': genomic_db, 'program':'blastn'}
    blast_runner = create_runner(tool='blast', parameters=parameters)
    blast_parser = get_alignment_parser('blast')
    filters = [{'kind'           : 'min_scores',
                'score_key'      : 'similarity',
                'min_score_value': 90,
               },
               {'kind'           : 'min_length',
                'min_length_bp'  : 20,
               }
              ]
    if not genomic_seqs_fpath:
        msg = 'No genomic sequence file defined for unique SNV filter'
        raise ValueError(msg)
    if not genomic_db:
        msg = 'No genomic blast database defined for unique SNV filter'
        raise ValueError(msg)
    genomic_seqs_fhand = open(genomic_seqs_fpath)
    genomic_seqs_index = SeqIO.index(genomic_seqs_fhand.name,
                                     guess_seq_file_format(genomic_seqs_fhand))

    def unique_contiguous_region_filter(sequence):
        '''It filters out the snv in regions repeated in the genome or
        discontiguous'''
        if sequence is None:
            return None

        for snv in sequence.get_features(kind='snv'):
            # Check if it is already done
            previous_result = _get_filter_result(snv, 'uniq_contiguous',
                                                 threshold=distance)
            if previous_result is not None:
                continue

            #we make a blast
            #with the sequence around the snv
            location = snv.location.start.position
            start = location - distance
            end = location + distance
            if start < 0:
                start = 0
            #print start, end
            seq_fragment = sequence[start:end]
            blast_fhand = blast_runner(seq_fragment)['blast']
            #now we parse the blast
            blast_result = blast_parser(blast_fhand)
            alignments = FilteredAlignmentResults(match_filters=filters,
                                              results=blast_result)
            #are there any similar sequences?
            try:
                alignment = alignments.next()
                result = True
            except StopIteration:
                #if there is no similar sequence we assume that is unique
                result = False
            if result:
                #how many matches, it should be only one
                num_hits = len(alignment['matches'])

                if num_hits > 1:
                    result = True
                else:
                    #how many match parts have the first match?
                    #we could do it with the blast result, but blast is not very
                    #good aligning, so we realign with est2genome
                    blast_fhand.seek(0)
                    sim_seqs = similar_sequences_for_blast(blast_fhand)
                    sim_seq = sim_seqs[0] if sim_seqs else None

                    introns = infer_introns_for_cdna(sequence=seq_fragment,
                                          genomic_seqs_index=genomic_seqs_index,
                                              similar_sequence=sim_seq,
                                              genomic_db=genomic_db)
                    if introns:
                        result = True
                    else:
                        result = False

            blast_fhand.close()
            _add_filter_result(snv, 'uniq_contiguous', result, distance)
        return sequence

    return unique_contiguous_region_filter

def create_close_to_intron_filter(distance):
    '''It returns a filter that filters snv by the proximity to introns.

    The introns should have been annotated before.'''

    def close_to_intron_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'close_to_intron',
                                                 threshold=distance)
            if previous_result is not None:
                continue

            location = snv.location.start.position
            result = False
            for intron in sequence.get_features(kind='intron'):
                if abs(location - intron.location.start.position) < distance:
                    result = True
            _add_filter_result(snv, 'close_to_intron', result,
                               threshold=distance)
        return sequence
    return close_to_intron_filter

def create_high_variable_region_filter(max_variability, window=0):
    'It returns a filter that filters snvs by region variability.'

    if window == 0:
        window = None

    def high_variable_region_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        snvs = list(sequence.get_features(kind='snv'))
        for snv in snvs:
            threshold = (max_variability, window)
            previous_result = _get_filter_result(snv, 'high_variable_reg',
                                                 threshold=threshold)
            if previous_result is not None:
                continue
            if window is None:
                snv_num = len(snvs)
                total_length = len(sequence)
            else:
                total_length = window
                snv_num = snvs_in_window(snv, snvs, window)
            variability = snv_num / float(total_length)
            if variability > max_variability:
                result = True
            else:
                result = False
            _add_filter_result(snv, 'high_variable_reg', result,
                               threshold=threshold)
        return sequence
    return high_variable_region_filter

def create_close_to_snv_filter(distance):
    '''It returns a filter that filters snv by the distance to other snvs.

    If the snv has another snv closer than DISTANCE, then this snv is
    filtered out'''
    def close_to_snv_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        snvs = list(sequence.get_features(kind='snv'))
        for snv in snvs:
            previous_result = _get_filter_result(snv, 'close_to_snv',
                                                 threshold=distance)
            if previous_result is not None:
                continue

            num_snvs = snvs_in_window(snv, snvs, distance * 2)
            if num_snvs > 1:
                result = True
            else:
                result = False

            _add_filter_result(snv, 'close_to_snv', result,
                               threshold=distance)
        return sequence
    return close_to_snv_filter

def create_snv_close_to_limit_filter(distance):
    '''This function is a function factory. This function is a filter than
     return true if those snvs are or not closer to the limit than max_distance
     '''
    def snv_close_to_reference_limit(sequence):
        '''True if the snv variation is close to the contig limits.

        It checks if the snv is not within the range covered by the
        consensus and if it's close to one of its limits. In both cases it will
        return True.
        '''
        if sequence is None:
            return None
        snvs = list(sequence.get_features(kind='snv'))
        for snv in snvs:
            previous_result = _get_filter_result(snv, 'close_to_limit',
                                                 threshold=distance)
            if previous_result is not None:
                continue
            location = int(snv.location.start.position)
            if location < distance or location + distance > len(sequence):
                result = True
            else:
                result = False
            _add_filter_result(snv, 'close_to_limit', result,
                               threshold=distance)
        return sequence
    return  snv_close_to_reference_limit

def create_major_allele_freq_filter(frequency, groups=None, group_kind=None):
    'It filters the snv in a seq by the frequency of the more frequent allele'
    def major_allele_freq_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'maf',
                                                 threshold=frequency)
            if previous_result is not None:
                continue
            maf = calculate_maf_frequency(snv, groups=groups,
                                          group_kind=group_kind)

            if maf > frequency or maf is None:
                result = True
            else:
                result = False
            _add_filter_result(snv, 'maf', result,
                               threshold=frequency)
        return sequence
    return major_allele_freq_filter

def create_kind_filter(kind):
    'It filters out the snvs with a different kind'
    def kind_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'by_kind', threshold=kind)
            if previous_result is not None:
                continue

            kind_ = calculate_snv_kind(snv)
            if kind != kind_:
                result = True
            else:
                result = False
            _add_filter_result(snv, 'by_kind', result, threshold=kind)
        return sequence
    return kind_filter

def create_cap_enzyme_filter(all_enzymes):
    'It filters the snv looking if it is detectable by rstriction enzymes'
    def cap_enzyme_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'cap_enzymes',
                                                 threshold=all_enzymes)
            if previous_result is not None:
                continue
            enzymes = calculate_cap_enzymes(snv, sequence,
                                            all_enzymes=all_enzymes)
            if len(enzymes) != 0:
                result = False
            else:
                result = True
            _add_filter_result(snv, 'cap_enzymes', result,
                               threshold=all_enzymes)
        return sequence
    return cap_enzyme_filter

def create_not_variable_in_group_filter(group_kind, groups, in_union=True):
    '''it filters looking if the list of reads is variable in the given
    conditions. It look in the'''

    if isinstance(groups, basestring):
        groups = (groups,)
    else:
        groups = tuple(groups)
    parameters = (group_kind, groups, in_union)

    def is_not_variable_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'is_not_variable',
                                                 threshold=parameters)
            if previous_result is not None:
                continue
            result = variable_in_groupping(group_kind, snv, groups, in_union,
                                           in_all_groups=False)
            result = bool(result)
            _add_filter_result(snv, 'is_not_variable', result,
                               threshold=parameters)
        return sequence

    return is_not_variable_filter

def create_is_variable_filter(group_kind, groups, in_union=True,
                              in_all_groups=True):
    '''it filters looking if the list of reads is variable in the given
    conditions. It look in the'''

    if isinstance(groups, basestring):
        groups = (groups,)
    else:
        groups = tuple(groups)
    parameters = (group_kind, groups, in_union)

    def is_variable_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'is_variable',
                                                 threshold=parameters)
            if previous_result is not None:
                continue
            result = variable_in_groupping(group_kind, snv, groups, in_union,
                                          in_all_groups=in_all_groups)
            result = True if result is None else not result

            _add_filter_result(snv, 'is_variable', result, threshold=parameters)
        return sequence

    return is_variable_filter

def create_min_groups_filter(min_groups, group_kind='read_groups'):
    'It filters snvs read in less groups (samples) than the min number given'

    parameters = (min_groups, group_kind)

    def min_groups_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'min_groups',
                                                 threshold=parameters)
            if previous_result is not None:
                continue

            #how many groups are in the alleles?
            groups = set()
            rg_info = snv.qualifiers['read_groups']
            for allele in snv.qualifiers['alleles'].values():
                read_groups = allele['read_groups'].keys()
                al_groups = [_get_group(rg, group_kind, rg_info) for rg in read_groups]
                groups = groups.union(set(al_groups))
            result = False if len(groups) >= min_groups else True

            _add_filter_result(snv, 'min_groups', result, threshold=parameters)
        return sequence

    return min_groups_filter
