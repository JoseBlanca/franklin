'''This module holds some utilities to build sequence cleaning pipelines.

In this module there is a collection of cleaning steps predefined. There
are steps that trim and mask the sequences. Every one of this step is a function
factory that will create the function that will do the actual job.

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

import logging, os, re
from itertools import tee
from Bio import SeqIO
import franklin
from franklin.utils.cmd_utils import create_runner
from franklin.seq.seqs import copy_seq_with_quality, Seq
from franklin.alignment_search_result import (FilteredAlignmentResults,
                                            get_alignment_parser)
from franklin.seq.seq_analysis import match_words

DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'data')

#The adaptors shorter than this length can not be processed by blast or
#exonerate, they should be processed by the word remover function
MIN_LONG_ADAPTOR_LENGTH = 20
TRIMMING_RECOMENDATIONS = 'trimming_recommendations'

def _add_trim_segments(segments, sequence, vector=True, trim=True):
    'It adds a segments or segments to the trimming recomendation in annotation'
    if not segments:
        return
    if 'trimming-recommendation' not in sequence.annotations:
        sequence.annotations[TRIMMING_RECOMENDATIONS] = {'vector':[],
                                                           'quality':[],
                                                           'mask':[]}

    if isinstance(segments, tuple):
        segments = [segments]
    trim_rec = sequence.annotations[TRIMMING_RECOMENDATIONS]
    if vector and trim:
        trim_rec['vector'].extend(segments)
    elif not vector and trim:
        trim_rec['quality'].extend(segments)
    elif not trim:
        trim_rec['mask'].extend(segments)

def _mask_sequence(sequence, segments):
    'It mask the given segments of the sequence'
    if not segments:
        return sequence
    segments = _get_all_segments(segments, len(sequence))
    seq = str(sequence.seq)
    new_seq = ''
    for segment in segments:
        start = segment[0][0]
        end   = 0 if segment[0][1] == 0 else  segment[0][1] + 1
        seq_  = seq[start:end]

        if segment[1]:
            seq_ = seq_.lower()
        new_seq += seq_
    return copy_seq_with_quality(sequence, seq=Seq(new_seq,
                                                   sequence.seq.alphabet))

def create_seq_trim_and_masker(mask=True, trim=True):
    'It actually trims the sequence taking into account trimming recomendations'
    def sequence_trimmer(sequence):
        'It trims the sequences'
        if not sequence:
            return None
        if not TRIMMING_RECOMENDATIONS in sequence.annotations:
            return sequence
        trim_rec = sequence.annotations[TRIMMING_RECOMENDATIONS]

        if mask:
            sequence = _mask_sequence(sequence, trim_rec.get('mask', []))
            if 'mask' in trim_rec:
                del(trim_rec['mask'])

        trim_locs = trim_rec.get('vector', []) + trim_rec.get('quality', [])
        if trim and trim_locs:
            trim_limits = _get_longest_non_matched_seq_region_limits(sequence,
                                                                     trim_locs)
            sequence = sequence[trim_limits[0]:trim_limits[1] + 1]
            if 'quality' in trim_rec:
                del(trim_rec['quality'])
            if 'vector' in trim_rec:
                del(trim_rec['vector'])
            if not mask and 'mask' in trim_rec:
                new_mask_segments = []
                for start, end in trim_rec['mask']:
                    start = start - trim_limits[0] if start - trim_limits[0] > 0 else 0
                    new_mask_segments.append((start, end - trim_limits[0]))
                trim_rec['mask'] = new_mask_segments
        sequence.annotations[TRIMMING_RECOMENDATIONS] = trim_rec
        return sequence

    return sequence_trimmer

def create_upper_mapper():
    'It returns a function that uppers the sequence'
    def upper(sequence):
        'Given a sequence it uppers it case'
        return sequence.upper()
    return upper

def create_edge_stripper(left_length=None, right_length=None):
    'It removes num of letters from seq.'
    def edge_stripper(sequence):
        'The real cleaner'
        if sequence is None:
            return None
        if left_length is None and right_length is None:
            return sequence
        segments = [(0, left_length - 1)] if left_length else []
        if right_length:
            seq_length = len(sequence)
            segments.append((seq_length - right_length  , seq_length - 1))

        _add_trim_segments(segments, sequence, vector=False)

        return sequence
    return edge_stripper

def create_striper_by_quality(quality_treshold, min_quality_bases=None,
                              min_seq_length=None, quality_window_width=None):
    '''It returns a function that given a sequence it returns a trimmed sequence
    without the bad quality region located at the extremes.

    This function calculates using an sliding window an average quality for
    every position. After that every position is classified, if its mean quality
    is above the given threshold is considered good, otherwise is bad. The
    regions to be trimmed are located looking for a contiguous number of bases
    with good quality.
    Arguments:
        - quality_threshold: Minimum quality to consider a position as good
        quality.
        - quality_window_width: The number of residues to use by the sliding
        window.
        - min_quality_bases: How many residues are necessary to define the
        start of a good quality zone. example:
                0 0 0 1 1 0 1 1 1
                      + +   > > >
                If the min_quality_bases is 2, the start end of the bad
                quality zone start with +. if it is 3, it starts with >
        - min_seq_length: Minimun sequence length of the trimmed sequence.
        If the length is less than that no sequence will be returned.
    '''
    def strip_seq_by_quality(sequence):
        '''Given a sequence it returns a trimmed sequence without the bad
        quality extremes.

        It could return None if the trimmed sequence is too short.
        '''
        if sequence is None:
            return None

        quality = sequence.qual
        seq_len = len(quality)
        # Define defaults
        min_quality_bases_, min_seq_length_, quality_window_width_ = \
                       _trim_seq_by_quality_defaults(seq_len, min_quality_bases,
                                                     min_seq_length,
                                                     quality_window_width)
        #calculate the sliding window qual
        qual_average_array = _calculate_sliding_window_qual(quality,
                                                          quality_window_width_)
        # Now we convert the quality average array in a boolean string. 1 if the
        # quality is good, and 0 if it is bad. Depending on the treshold
        boolean_quality_treshold = _quality_to_boolean(quality_treshold,
                                                       qual_average_array)

        start, end = _trim_bad_qual_extremes(boolean_quality_treshold,
                                             min_quality_bases_)

        segments = [(0, start -1), (end + 1, len(sequence) -1)]
        _add_trim_segments(segments, sequence, vector=False)
        #print 'seq', sequence.seq
        if (end - start) < min_seq_length_:
            return None
        else:
            return sequence
    return strip_seq_by_quality

def _trim_bad_qual_extremes(bool_seq, min_quality_bases):
    '''It returns start and and of the new sequence. Givig the 0/1 string.'''

    def _get_good_qual_extreme_loc(bool_seq, min_quality_bases, extreme=None):
        'It look qood quality extreme start. Good quality is min_qual_bases'
        extreme_pos     = 0
        good_count      = 0
        in_good_section = True
        if extreme == 'end':
            bool_seq = bool_seq[::-1]

        for i, bool_ in enumerate(bool_seq):
            if bool_ == '1' and in_good_section:
                good_count += 1
            elif bool_ == '1' and not in_good_section:
                extreme_pos = i
                in_good_section = True
            elif bool_ == '0':
                in_good_section = False
            if good_count == min_quality_bases:
                break

        if extreme == 'end':
            extreme_pos = len(bool_seq) - extreme_pos

        return extreme_pos

    start = _get_good_qual_extreme_loc(bool_seq, min_quality_bases, 'start')
    end   = _get_good_qual_extreme_loc(bool_seq, min_quality_bases, 'end') - 1
    return start, end

def _calculate_sliding_window_qual(quality, quality_window_width):
    '''It takes into account the neighbour qualities to calculate the new
    quality. The window size is variable '''
    average_quality = []
    #pylint: disable-msg=W0612
    for i, single_qual in enumerate(quality):
        if quality_window_width > i:
            qual_slice = quality[0: i + quality_window_width +1]
        else:
            qual_slice = quality[i -quality_window_width:
                                 i + quality_window_width +1]
        average = sum(qual_slice) / float(len(qual_slice))
        average_quality.append(average)
    return average_quality

def _quality_to_boolean(quality_treshold, qual_average_array):
    '''It converts the average_quality array in a boolean string.  1 if the
    quality is good, and 0 if it is bad. Depending on the treshold'''
    boolean_qual = []
    for quality in qual_average_array:
        if quality > quality_treshold:
            boolean_qual.append('1')
        else:
            boolean_qual.append('0')
    return ''.join(boolean_qual)

def _trim_seq_by_quality_defaults(seq_len, min_quality_bases, min_seq_length,
                                  quality_window_width):
    '''We define defaults for the  trim_seq_by_quality function.

    Using comments of Giles Weaver <giles.weaver@googlemail.com> to the
    biopython list.

     '''
    if min_quality_bases is None:
        if seq_len < 30:
            min_quality_bases = 3
        elif seq_len < 45:
            min_quality_bases = 3
        elif seq_len < 70:
            min_quality_bases = 6
        elif seq_len < 300:
            min_quality_bases = 15
        else:
            min_quality_bases = 25

    if quality_window_width is None:
        if seq_len < 30:
            quality_window_width = 1
        if seq_len < 70:
            quality_window_width = 2
        elif seq_len < 300:
            quality_window_width = 3
        else:
            quality_window_width = 4

    if min_seq_length is None:
        fiftypercent = seq_len / 2
        if fiftypercent > 100:
            min_seq_length = 100
        else:
            min_seq_length = fiftypercent
    return (min_quality_bases, min_seq_length, quality_window_width)

def create_masker_for_low_complexity():
    'It creates a masker function for low complexity sections that uses mdust'
    mask_low_complex_by_seq = create_runner(tool='mdust')

    def mask_low_complexity(sequence):
        '''It adds a mask to the sequence where low complexity is found

        It uses mdust from the seqclean package
        '''
        if sequence is None:
            return None
        fhand = mask_low_complex_by_seq(sequence)['sequence']
        segments = []
        for line in fhand:
            start, end = line.strip().split()[-2:]
            segments.append((int(start) - 1, int(end) - 1))
        _add_trim_segments(segments, sequence, trim=False)

        return sequence
    return mask_low_complexity

def create_masker_for_polia():
    'It creates a masker function that will mask poly-A tracks'
    parameters = {'min_score':'10', 'end':'x', 'incremental_dist':'20',
                      'fixed_dist':None}
    mask_polya_by_seq = create_runner(tool='trimpoly', parameters=parameters)
    def mask_polya(sequence):
        '''It adds a mask to the sequence where the poly-A is found.

        It uses trimpoly from seqclean package
        '''
        if sequence is None:
            return None

        fhand = mask_polya_by_seq(sequence)['sequence']
        segments = _segments_from_trimpoly(fhand, sequence)
        _add_trim_segments(segments, sequence, trim=False)
        return sequence
    return mask_polya

def create_word_masker(words, beginning=True):
    'It removes the given words if they are in the start of the seq'
    if not beginning:
        words = [re.compile(word) for word in words]

    def word_remover(sequence):
        'The remover'
        if sequence is None:
            return None
        str_seq = str(sequence.seq)
        for word in words:
            segment = None
            if beginning and str_seq.startswith(word):
                segment = (0, len(word) - 1)
            elif not beginning:
                segment = _findallmatches(word, str_seq)

            if segment:
                _add_trim_segments(segment, sequence, trim=False)
        return sequence

    return word_remover

def _findallmatches(regex, string):
    '''Given a compiled regex it finds all matches in the string and return
    their location'''
    matches = regex.finditer(string)
    if matches:
        return [(match.start(), match.end() - 1) for match in matches]
    else:
        return []

def create_striper_by_quality_trimpoly(ntrim_above_percent=2):
    '''It creates a function that removes bad quality regions.

    It uses trimpoly's quality checks.'''
    ntrim_above_percent = '%.1f' % ntrim_above_percent

    def strip_seq_by_quality_trimpoly(sequence):
        '''It strips the sequence where low quality is found

        It uses trimpoly from seqclean package.
        This program does not work well with short sequences.
        '''
        if sequence is None:
            return None

        if len(sequence) < 80:
            return None

        parameters = {'only_n_trim':None,
                      'ntrim_above_percent': ntrim_above_percent}
        mask_polya_by_seq = create_runner(tool='trimpoly',
                                          parameters=parameters)
        fhand = mask_polya_by_seq(sequence)['sequence']
        segments = _segments_from_trimpoly(fhand, sequence)
        _add_trim_segments(segments, sequence, vector=False)
        return sequence
    return strip_seq_by_quality_trimpoly

def _segments_from_trimpoly(fhand_trimpoly_out, sequence):
    '''It return new sequence giving that trimpoly output and the old sequence
     trim option is used to trim os mask the low quality sequence '''

    trimp_data = fhand_trimpoly_out.readline().split()
    end5 = int(trimp_data[2]) - 1
    end3 = int(trimp_data[3])
    segments = [(0, end5 - 1)] if end5 else []
    if end3 < len(sequence):
        segments.append((end3, len(sequence)))

    return segments

def _lucy_mapper(sequence, index):
    '''It processes the sequence taking the lucy result.

    The index is the lucy result indexed
    '''
    name = sequence.name
    lucy_seq = index[name]
    lucy_description = lucy_seq.description

    if lucy_description is None:
        return sequence

    start, end = lucy_description.split()[-2:]

    #Now we have the segment in 0 coord system and into a segment
    matched_segment = (int(start) - 1, int(end) -1)

    # we need the non matched segment
    segments = _get_non_matched_from_matched_locations([matched_segment],
                                                       len(sequence))

    # add this to the sequence annotation
    _add_trim_segments(segments, sequence)
    return sequence

def create_striper_by_quality_lucy(parameters=None):
    '''It creates a function that removes bad quality regions using lucy.

    The function will take a sequence iterator and it will return a new sequence
    iterator with the processed sequences in it.'''
    run_lucy_for_seqs = create_runner(tool='lucy', parameters=parameters)

    def strip_seq_by_quality_lucy(sequences):
        '''It trims the bad quality regions from the given sequences.

        It uses lucy external program. It returns a sequence iterator and a,
        list of output files. (A seq file and a qual file). We return the files
        because the iterator feeds from them and they are temporary files that
        will be removed as soon as they get out of scope.
        '''
        #pylint: disable-msg=W0612
        #now we run lucy
        sequences, sequences_for_lucy = tee(sequences, 2)
        seq_out_fhand = run_lucy_for_seqs(sequences_for_lucy)['sequence'][0]

        # index the lucy result
        result_index = SeqIO.index(seq_out_fhand.name, "fasta")

        # process each sequence and
        for sequence in sequences:
            yield _lucy_mapper(sequence, result_index)

    return strip_seq_by_quality_lucy

#pylint:disable-msg=C0103
def create_vector_striper_by_alignment(vectors, aligner):
    '''It creates a function which will remove vectors from the given sequence.

    It looks for the vectors comparing the sequence with a vector database. To
    do these alignments two programs can be used, exonerate and blast. Exonerate
    requires a fasta file with the vectors and blast and indexed blast database.
    '''
    #depending on the aligner program we need different parameters and filters
    parameters = {'exonerate': {'target':vectors},
                  'blast'    : {'database': vectors, 'program':'blastn'},
                  'blast+'    : {'database': vectors, 'program':'blastn'}}

    filters = {'exonerate': [{'kind'          : 'min_scores',
                             'score_key'      : 'similarity',
                             'min_score_value': 96},
                             {'kind'          : 'min_length',
                              'min_length_bp' : MIN_LONG_ADAPTOR_LENGTH}],
               'blast':      [{'kind'         : 'min_scores',
                             'score_key'      : 'similarity',
                             'min_score_value': 96},
                             {'kind'          : 'min_length',
                              'min_length_bp' : MIN_LONG_ADAPTOR_LENGTH}],
               'blast+':      [{'kind'         : 'min_scores',
                             'score_key'      : 'similarity',
                             'min_score_value': 96},
                             {'kind'          : 'min_length',
                              'min_length_bp' : MIN_LONG_ADAPTOR_LENGTH}]}

    aligner_ = create_runner(tool=aligner, parameters=parameters[aligner])
    parser   = get_alignment_parser(aligner)

    def strip_vector_by_alignment(sequence):
        '''It strips the vector from a sequence.

        It returns a striped sequence with the longest segment without vector.
        '''

        if sequence is None:
            return None
        if vectors is None:
            return sequence

        # first we are going to align he sequence with the vectors
        alignment_fhand = aligner_(sequence)[aligner]
        # We need to parse the result
        alignment_result = parser(alignment_fhand)

        # We filter the results with appropriate  filters
        alignments = FilteredAlignmentResults(match_filters=filters[aligner],
                                              results=alignment_result)

        alignment_matches = _get_non_matched_locations(alignments)
        alignment_fhand.close()
        segments  = _get_longest_non_matched_seq_region_limits(sequence,
                                                              alignment_matches)
        if segments is None:
            return None

        segments  = _get_non_matched_from_matched_locations([segments],
                                                            len(sequence))
        _add_trim_segments(segments, sequence)
        return sequence

    return strip_vector_by_alignment

def create_word_striper_by_alignment(words):
    '''It creates a function which will remove words from the given sequence.

    It matches the words against the sequence and leaves the longest non-matched
    part (unless is the first part of the sequence).
    '''

    def strip_words_by_matching(sequence):
        '''It strips the given words from a sequence.

        It returns a striped sequence with the longest segment without the
        words.
        '''
        if sequence is None:
            return None
        if not words:
            return sequence

        alignments = match_words(sequence, words)
        if not alignments:
            return sequence
        locations = _get_non_matched_locations(alignments)
        segments  = _get_longest_non_matched_seq_region_limits(sequence,
                                                              locations)
        if segments is None:
            return None
        segments  = _get_non_matched_from_matched_locations([segments],
                                                           len(sequence))
        _add_trim_segments(segments, sequence)
        return sequence

    return strip_words_by_matching


def _merge_overlaping_locations(locations):
    '''It merges overlaping locations

    It accept a list of (start, end) tuples.
    '''
    #we collect all start and ends in a list marking if they are start or end
    limits = [] #all starts and ends
    for location in locations:
        limit_1 = {
            'type' : 'start',
            'location' : location[0]
        }
        limit_2 = {
            'type' : 'end',
            'location' : location[1]
        }
        limits.append(limit_1)
        limits.append(limit_2)

    def cmp_limits(limit1, limit2):
        'It returns -1, 0 or 1 depending on the limit order'
        loc1 = limit1['location']
        loc2 = limit2['location']
        #we sort using its locations
        if loc1 != loc2:
            return loc1 -loc2
        #unless they're in the same location.
        type1 = limit1['type']
        type2 = limit2['type']
        if type1 == type2:
            return 0
        elif type1 == 'start' and type2 == 'end':
            return -1
        else:
            return 1

    #we sort the start and ends by the start position
    limits = sorted(limits, cmp_limits)

    #if there are contiguous locations we should merged them
    new_limits = []
    merge = False
    for index, limit in enumerate(limits):
        try:
            next_limit = limits[index + 1]
        except IndexError:
            next_limit = None
        if (next_limit and
            limit['type'] == 'end' and
            next_limit['type'] == 'start' and
            limit['location'] + 1 == next_limit['location']):
            merge = True
            continue
        if merge:
            merge = False
            continue
        new_limits.append(limit)
    limits = new_limits

    #now we create the merged locations
    starts = 0
    merged_locations = []
    for limit in limits:
        if limit['type'] == 'start':
            starts += 1
            if starts == 1:
                start = limit['location']
        elif limit['type'] == 'end':
            starts -= 1
            if starts == 0:
                end = limit['location']
                merged_locations.append((start, end))

    return merged_locations

def _get_non_matched_locations(alignments):
    '''It detects query's non matched fragments in an alignment

    It returns a list of  (start, end) tuples
     '''
    # we need the matches sorted by its start and merged when two overlap
    # we collect all match starts and ends in tuples
    locations = []
    for alignment in alignments:
        for match in alignment['matches']:
            for match_part in match['match_parts']:
                locations.append((match_part['query_start'],
                                  match_part['query_end']))
    return locations

def _get_unmasked_locations(seq):
    '''It detects the unmasked regions of a sequence

    It returns a list of (start, end) tuples'''

    locations = []
    in_masked_section = True
    for i, letter in enumerate(seq):
        unmasked =  str(letter).isupper()
        if unmasked and in_masked_section:
            start = i
            in_masked_section = False
        elif not unmasked and not in_masked_section:
            end = i - 1
            locations.append((start, end))
            in_masked_section = True
    else:
        if unmasked and not in_masked_section:
            end = i
            locations.append((start, end))

    return locations

def _get_all_segments(segments, seq_len):
    '''Given a set of some non overlaping regions it returns all regions.

    input:  --- ----       ----
    output: ---+----+++++++----

    '''
    if not segments:
        return [((0, seq_len - 1), False)]

    all_segments = []
    #If the first segment doesn't start at zero we create a new one
    if segments[0][0] == 0:
        start = segments[0][1] + 1
        all_segments.append((segments.pop(0), True))
    else:
        start = 0

    for loc in segments:
        all_segments.append(((start, loc[0] - 1), False))
        all_segments.append((loc, True))
        start = loc[1] + 1
    else:
        #if the last segment does not ends at the end of the sequence we add
        #an extra one
        end = seq_len - 1
        if start <= end:
            all_segments.append(((start, end), False))
    return all_segments

def _get_non_matched_from_matched_locations(segments, seq_len):
    'Given a set of regions in a seq it returns the complementary ones'
    non_matched = []
    for segment in _get_all_segments(segments, seq_len):
        if not segment[1]:
            non_matched.append(segment[0])
    return non_matched


def _get_longest_non_matched_seq_region_limits(seq, locations):
    'Given a seq and the locations it returns the longest region limits'
    #they are one based, we do them 0 based
    #locations = [(loc[0] - 1, loc[1] - 1) for loc in locations]

    # Once we collect all the start and ends in a list of tuples, we are going
    # to merge overlapings
    locations = _merge_overlaping_locations(locations)

    #we want the non-matched locations
    locations = _get_non_matched_from_matched_locations(locations, len(seq))

    #if there is one that starts at zero we remove it, because it should be
    #vector
    #if locations[0][0] == 0:
    #    locations.pop(0)

    #now we return the longest one
    longest_loc = None
    for loc in locations:
        if not longest_loc:
            longest_loc = loc
            continue
        if longest_loc[1] - longest_loc[0] < loc[1] - loc[0]:
            longest_loc = loc
    return longest_loc

def _get_matched_locations(seq, locations, min_length):
    'It returns a seq iterator from a seq. To split the seq it uses the'

    for i, (start, end) in enumerate(locations):
        seq1 = seq.seq[start:end+1]
        if len(seq1) < min_length:
            continue
        if seq.qual is not None:
            qual = seq.qual[start:end+1]
        else:
            qual = None
        if seq.name is not None:
            name = '%s_%d' % (seq.name, i + 1)
        yield copy_seq_with_quality(seq, seq=seq1, qual=qual, name=name)

def split_seq_by_masked_regions(seq_iter, min_length):
    '''It takes a sequence iterator and return another iterator of sequences.
    The masked section of the sequences have been removed from these sequences
    and the resulting seqs are returned as new seqs.
        example:  AATTAATTGGagatagatAATTGATGAATGAtaaaaaaaaGATAGATAGAGAGT
            newseq1 = AATTAATTGG
            newseq2 = AATTGATGAATGA
            newseq3 = GATAGATAGAGAGT
    '''
    for seq in seq_iter:
        locations = _get_unmasked_locations(seq)
        # if there are no masked sections, we return the whole sequence
        if not locations:
            yield seq
            continue
        seqs = _get_matched_locations(seq, locations, min_length)
        for new_seq in seqs:
            yield new_seq


