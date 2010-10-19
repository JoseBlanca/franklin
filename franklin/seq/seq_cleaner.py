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

import logging, os
from tempfile import NamedTemporaryFile
from itertools import tee

import franklin
from franklin.utils.cmd_utils import (create_runner,
                                      run_repeatmasker_for_sequence)
from franklin.seq.writers import fasta_str
from franklin.seq.readers import seqs_in_file
from franklin.utils.seqio_utils import get_content_from_fasta
from franklin.seq.seqs import copy_seq_with_quality, Seq
from franklin.alignment_search_result import (FilteredAlignmentResults,
                                            get_alignment_parser)
from franklin.seq.seq_analysis import match_words

DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'data')

#The adaptors shorter than this length can not be processed by blast or
#exonerate, they should be processed by the word remover function
MIN_LONG_ADAPTOR_LENGTH = 15

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
        rigth_limit = len(sequence) - right_length
        return sequence[left_length:rigth_limit]
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
        new_seq = sequence[start:end]
        if len(new_seq.qual) < min_seq_length_:
            return None
        else:
            return new_seq
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
    end   = _get_good_qual_extreme_loc(bool_seq, min_quality_bases, 'end')
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
        #pylint:disable-msg=W0612
        seq = get_content_from_fasta(fhand)[-1]
        return copy_seq_with_quality(sequence, seq=Seq(seq,
                                                       sequence.seq.alphabet))
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
        return _sequence_from_trimpoly(fhand, sequence, trim=False)
    return mask_polya

def create_word_masker(words, beginning=True):
    'It removes the given words if they are in the start of the seq'
    def word_remover(sequence):
        'The remover'
        if sequence is None:
            return None
        str_seq = str(sequence.seq)
        for word in words:
            if beginning and str_seq.startswith(word):
                seq = str_seq[:len(word)].lower() + str_seq[len(word):]
                seq = Seq(seq)
                sequence = copy_seq_with_quality(sequence, seq=seq)
            elif not beginning:
                seq = str_seq.replace(word, word.lower())
                seq = Seq(seq)
                sequence = copy_seq_with_quality(sequence, seq=seq)

        return sequence

    return word_remover

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
        return _sequence_from_trimpoly(fhand, sequence, trim=True)
    return strip_seq_by_quality_trimpoly

def _sequence_from_trimpoly(fhand_trimpoly_out, sequence, trim):
    '''It return new sequence giving that trimpoly output and the old sequence
     trim option is used to trim os mask the low quality sequence '''

    trimp_data = fhand_trimpoly_out.readline().split()
    end5 = int(trimp_data[2]) - 1
    end3 = int(trimp_data[3])
    new_sequence = ''
    str_seq = str(sequence.seq)
    if not trim:
        new_sequence  += str_seq[:end5].lower()
    new_sequence += str_seq[end5:end3]

    if trim:
        return sequence[end5:end3]
    else:
        new_sequence += str_seq[end3:].lower()
        return copy_seq_with_quality(sequence, seq=Seq(new_sequence,
                                                       sequence.seq.alphabet))

def create_striper_by_quality_lucy():
    'It creates a function that removes bad quality regions using lucy'
    #This function has been deprecated because when running using condor lucy
    #tends to get stuck without apparent reason. We think that it might be due
    #to the fact that lucy was created to analyze big fasta files with a lot
    #of sequences in them. In this case we were using lucy to analyze individual
    #sequences.
    #To workaround this problem we have created a new function that analyze
    #all the sequences at ones.
    msg  = "DEPECRATED, Lucy might get stuck when using this function."
    msg += "Use create_striper_by_quality_lucy2"
    print msg
    def strip_seq_by_quality_lucy(sequence):
        '''It trims from the sequence  bad quality sections.

        It uses lucy external program
        '''
        if sequence is None:
            return None

        if sequence.name is None:
            raise ValueError('lucy requires that the sequence has a name')
        elif len(sequence) < 130:
            if sequence.name is None:
                print_name = ''
            else:
                print_name = sequence.name
            msg = 'lucy: %s is shorter than 130. lucy min seq length 130' % \
                                                                    print_name
            logging.warning(msg)
            return None
        #we run lucy
        run_lucy_for_seq = create_runner(tool='lucy')
        #pylint: disable-msg=W0612
        seq_out_fhand = run_lucy_for_seq(sequence)[0]

        #from the output we know where to strip the seq
        description = get_content_from_fasta(seq_out_fhand)[1]

        #lucy can consider that the seq is low qual and return an empty file
        if description is None:
            return None
        start, end = description.split()[-2:]
        #we count from zero
        start, end = int(start) - 1, int(end)
        striped_seq = sequence[start:end]
        return striped_seq
    return strip_seq_by_quality_lucy

def _extract_name_description(seq_iter):
    '''using a sequence iterator, it creates a dictionary with name as key an
    description as value.'''
    name_description = {}
    for seq in seq_iter:
        try:
            name        = seq.name
            description = seq.description
        except AttributeError:
            pass
        name_description[name] = description

    return name_description

def create_striper_by_quality_lucy2(parameters=None):
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
        sequences, sequences_copy = tee(sequences, 2)
        seq_out_fhand, qual_out_fhand = run_lucy_for_seqs(sequences)['sequence']
        # we need to preserve the original description field.
        descriptions = _extract_name_description(sequences_copy)

        #now we have to clean all sequences with the description found on the
        #output seq file
        seq_iter = seqs_in_file(seq_out_fhand, qual_out_fhand)
        striped_seq_fhand = NamedTemporaryFile(suffix='.seq.fasta')
        striped_qual_fhand = NamedTemporaryFile(suffix='.qual.fasta')
        for seq in seq_iter:
            description = seq.description
            #lucy can consider that the seq is low qual and return an empty file
            if description is None:
                striped_seq = None
            start, end = description.split()[-2:]
            #we count from zero
            start, end = int(start) - 1, int(end)
            striped_seq = seq[start:end]
            # recover original description
            desc_orig = descriptions[seq.name]
            seq_seq = striped_seq.seq
            seq_qual = striped_seq.qual
            if len(seq_seq) != len(seq_qual):
                msg = 'Malformed lucy output, name: %s, seq: %s, qual:%s' % \
                                 (striped_seq.name, str(seq_seq), str(seq_qual))
                raise RuntimeError(msg)
            seq_str = fasta_str(seq_seq, striped_seq.name, desc_orig)
            striped_seq_fhand.write(seq_str)
            stripped_qual = ' '.join([str(q_val) for q_val in seq_qual])
            qual_str = fasta_str(stripped_qual, striped_seq.name, desc_orig)
            striped_qual_fhand.write(qual_str)
        seq_out_fhand.close()
        qual_out_fhand.close()
        #os.remove(seq_out_fhand.name)
        #os.remove(qual_out_fhand.name)
        striped_seq_fhand.flush()
        striped_qual_fhand.flush()
        seq_iter = seqs_in_file(striped_seq_fhand, striped_qual_fhand)
        return seq_iter, [striped_seq_fhand, striped_qual_fhand]
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
        return _get_longest_non_matched_seq_region(sequence, alignment_matches)

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
        return _get_longest_non_matched_seq_region(sequence, locations)

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

def _get_non_matched_from_matched_locations(locations, seq_len):
    'Given a set of regions in a seq it returns the complementary ones'
    if not locations:
        return [(0, seq_len - 1)]

    comp_locations = []
    #the first complementary region starts at 0 unless the first given location
    #starts at 0
    if locations[0][0] == 0:
        start = locations[0][1] + 1
        locations.pop(0)
    else:
        start = 0
    for loc in locations:
        comp_locations.append((start, loc[0] - 1))
        start = loc[1] + 1
    else:
        #if the last location do not end at the end of the seq we add the
        #last complementary location
        end = seq_len - 1
        if start <= end:
            comp_locations.append((start, end))
    return comp_locations

def _get_longest_non_matched_seq_region(seq, locations):
    'Given a seq and the locations it returns the longest region'

    # If no locations we return the complete sequence location
    if not locations:
        return seq

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

    if longest_loc:
        return seq[longest_loc[0]:longest_loc[1] + 1]
    else:
        return None

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

def create_masker_repeats_by_repeatmasker(species='eudicotyledons'):
    '''It creates a function that mask repeats from a sequence'''

    def mask_repeats_by_repeatmasker(sequence):
        '''It returns a sequence with the repetitive and low complex elements
         masked'''
        if sequence is None:
            return None

        seq_fhand = run_repeatmasker_for_sequence(sequence, species)
        #pylint: disable-msg=W0612
        masked_seq = get_content_from_fasta(seq_fhand)[2]
        return copy_seq_with_quality(sequence, seq=Seq(masked_seq,
                                                       sequence.seq.alphabet))
    return mask_repeats_by_repeatmasker

def create_short_adaptor_striper(adaptors):
    '''It creates a function which removes the adaptors from the given sequence.

    The adaptors should be a list. It will look for them using regular
    expressions.
    '''
    def strip_vector_by_alignment(sequence):
        '''It strips the vector from a sequence.

        It returns a striped sequence with the longest segment without vector.
        '''
        if sequence is None:
            return None
        if adaptors is None:
            return sequence
        alignments = parser(alignment_fhand)
        alignment_matches = _get_non_matched_locations(alignments)
        return _get_longest_non_matched_seq_region(sequence, alignment_matches)

    return strip_vector_by_alignment

