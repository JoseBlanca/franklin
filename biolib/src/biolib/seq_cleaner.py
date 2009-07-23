'''
Functions to clean sequences.

All functions here have the same interface. You pass a SeqWithQuality to the
fucntion and it returns another SeqWithQuality.
Depending on the algorith used it mask or it trim the bad seq/quality section.

As bad we undertand low quality section, low complexity section, poliA section,
...
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
from itertools import imap, ifilter, tee
import re, logging
from tempfile import NamedTemporaryFile

from biolib.biolib_cmd_utils import create_runner, run_repeatmasker_for_sequence
from biolib.biolib_utils import (get_content_from_fasta, seqs_in_file,
                                 get_safe_fname, fasta_str)
from biolib.seqs import  SeqWithQuality
from biolib.seq_filters import create_length_filter
from biolib.alignment_search_result import  (FilteredAlignmentResults,
                                             get_alignment_parser)



def create_striper_by_quality(quality_treshold, min_quality_bases=None,
                              min_seq_length=None, quality_window_width=None):
    'It creates a stripper that trims from the sequence  bad quality extremes.'
    def strip_seq_by_quality(sequence):
        '''It trims from the sequence  bad quality extremes.

        This function uses our algorithm based in the average of the neighbours
        qualities. It calculates the average quality of each position taking
        into account the qualities of the neighbours.
        Arguments:
            .- quality_threshold. Minimun quality needed to pass the check. For
            phred qualities a number of 20 is a good threshold
            .- quality_window_width: The number of nucleotide qualities of each
             side to calculate the average.
            .- min_quality_bases: How many nucleotides we need to define the
            start of a good quality zone. example:
                    0 0 0 1 1 0 1 1 1
                          + +   > > >
                    I the min_quality_bases is 2, the start end of the bad
                    quality zone start with +. if it is 3, it starts with >
            .- min_seq_length: Minimun seq length of the new sequence to mark as
            valid. Otherwise we dismiss the sequence
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
    start = re.search('0*(1{%d})' % min_quality_bases, bool_seq).start(1)

    #First I reverse de quality string
    bool_seq = bool_seq[::-1]
    end = re.match('0*(1{%d})' % min_quality_bases , bool_seq).start(1)
    end = len(bool_seq) - end
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
        if seq_len < 70:
            min_quality_bases = 6
        elif seq_len < 300:
            min_quality_bases = 15
        else:
            min_quality_bases = 25

    if quality_window_width is None:
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
    'It creates a masker function for low complexity sections'
    mask_low_complex_by_seq = create_runner(kind='mdust')

    def mask_low_complexity(sequence):
        '''It adds a mask to the sequence where low clomplexity is found

        It uses mdust from the seqclean package
        '''
        if sequence is None:
            return None
        fhand = mask_low_complex_by_seq(sequence)[0]
        #pylint:disable-msg=W0612
        name, desc, seq = get_content_from_fasta(fhand)
        return sequence.copy(seq=seq)
    return mask_low_complexity

def create_masker_for_polia():
    'It creates a masker function for polia slices'
    parameters = {'min_score':'10', 'end':'x', 'incremental_dist':'20',
                      'fixed_dist':None}
    mask_polya_by_seq = create_runner(kind='trimpoly', parameters=parameters)
    def mask_polya(sequence):
        '''It adds a mask to the sequence where poly A are found

        It uses trimpoly from seqclean package
        '''
        if sequence is None:
            return None

        fhand = mask_polya_by_seq(sequence)[0]
        return _sequence_from_trimpoly(fhand, sequence, trim=False)
    return mask_polya

def create_striper_by_quality_trimpoly():
    'It creates a function hat stripes seq using trimpoly\'s quality checks '


    def strip_seq_by_quality_trimpoly(sequence):
        '''It strip the sequence where low quality is found

        It uses trimpoly from seqclean package.
        This program does not work well with short sequences, and it only can
        look at the last 30 nucleotides max
        '''
        if sequence is None:
            return None

        if len(sequence) < 80:
            if sequence.name is None:
                print_name = ''
            else:
                print_name = sequence.name
            msg = 'trimpoly: Sequence %s shorter than 80 nt removed' % \
                                                                    print_name
            logging.warning(msg)
            return None

        parameters = {'only_n_trim':None, 'ntrim_above_percent': '3'}
        mask_polya_by_seq = create_runner(kind='trimpoly',
                                          parameters=parameters)
        fhand = mask_polya_by_seq(SeqWithQuality(seq=sequence))[0]
        return _sequence_from_trimpoly(fhand, sequence, trim=True)
    return strip_seq_by_quality_trimpoly

def _sequence_from_trimpoly(fhand_trimpoly_out, sequence, trim):
    '''It return new sequence giving tha trimpoly output and the old sequence
     trim option is used to trim os mask the low quality sequence '''

    trimp_data = fhand_trimpoly_out.readline().split()
    end3 = int(trimp_data[2]) - 1
    end5 = int(trimp_data[3])
    new_sequence = ''
    str_seq = sequence.seq
    if not trim:
        new_sequence  += str_seq[:end3].lower()
    new_sequence += str_seq[end3:end5]
    if not trim:
        new_sequence += str_seq[end5:].lower()
    if trim:
        return sequence[end3:end5]
    else:
        seq_class = sequence.seq.__class__
        return sequence.copy(seq=seq_class(new_sequence))

def create_striper_by_quality_lucy():
    'It creates a function hat stripes seq using lucys quality checks '
    msg  = "DEPECRATED, This function is deprecated. "
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
        run_lucy_for_seq = create_runner(kind='lucy')
        #pylint: disable-msg=W0612
        seq_out_fhand, qual_out_fhand = run_lucy_for_seq(sequence)

        #from the output we know where to strip the seq
        name, description, seq = get_content_from_fasta(seq_out_fhand)

        #lucy can consider that the seq is low qual and return an empty file
        if description is None:
            return None
        start, end = description.split()[-2:]
        #we count from zero
        start, end = int(start) - 1, int(end)
        striped_seq = sequence[start:end]
        return striped_seq
    return strip_seq_by_quality_lucy

def create_striper_by_quality_lucy2():
    'It creates a function hat stripes seq using lucys quality checks '
    def strip_seq_by_quality_lucy(sequences):
        '''It trims from the sequences  bad quality sections.

        It uses lucy external program. It return a sequence iterator and a,
        list of output files. (A seq file and a qual file). We return the files
        because the iterator use them and as they are temporal files, we need
        to return them in order to be usefull for the iterator
        '''
        #we prepare the function that will run lucy
        run_lucy_for_seqs = create_runner(kind='lucy')
        #pylint: disable-msg=W0612
        #now we run lucy
        seq_out_fhand, qual_out_fhand = run_lucy_for_seqs(sequences)

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
            seq_str = fasta_str(striped_seq.seq, striped_seq.name)
            striped_seq_fhand.write(seq_str)
            stripped_qual = ' '.join([str(q_val) for q_val in striped_seq.qual])
            qual_str = fasta_str(stripped_qual, striped_seq.name)
            striped_qual_fhand.write(qual_str)
        seq_iter = seqs_in_file(striped_seq_fhand, striped_qual_fhand)
        return seq_iter, [striped_seq_fhand, striped_qual_fhand]
    return strip_seq_by_quality_lucy


#pylint:disable-msg=C0103
def create_vector_striper_by_alignment(vectors, aligner):
    'It creates a function witch we can pass a sequence to search from vectors'
    #depending on the aligner program we need different parameters and filters
    parameters = {'exonerate': {'target':vectors},
                  'blast'    : {'database': vectors, 'program':'blastn'}}

    filters = {'exonerate': [{'kind'          : 'min_scores',
                             'score_key'      : 'similarity',
                             'min_score_value': 96},
                             {'kind'          : 'min_length',
                              'min_length_bp' : 15}],

               'blast':      [{'kind'         : 'min_scores',
                             'score_key'      : 'similarity',
                             'min_score_value': 96},
                             {'kind'          : 'min_length',
                              'min_length_bp' : 15}]}
    aligner_ = create_runner(kind=aligner, parameters=parameters[aligner])
    parser   = get_alignment_parser(aligner)

    def strip_vector_by_alignment(sequence):
        '''It strips the vector from a sequence

        It returns a striped sequence with the longest segment without vector.
        It can work with  blast or exonerate. And takes vector information from
        a database or a fasta file
        '''
        if sequence is None:
            return None
        #
        # first we are going to align he sequence with the vectors
        alignment_fhand = aligner_(sequence)[0]

        # We need to parse the result
        alignment_result = parser(alignment_fhand)

        # We filter the results with appropiate  filters
        alignments = FilteredAlignmentResults(match_filters=filters[aligner],
                                              results=alignment_result)

        start_end = _get_longest_non_matched_fragment(alignments, sequence)
        if start_end is None:
            return None
        else:
            start, end = start_end
            return sequence[start:end]
    return strip_vector_by_alignment

def _merge_overlaping_locations(locations):
    '''It merges overlaping locations

    It accept a list of (start, end) tuplas.
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

#pylint:disable-msg=C0103
def _get_longest_non_matched_fragment(alignments, query):
    '''It get the longest non matched fragment from matches result

    It use salignment result dict and after merging overlaping matches, it looks
    for the longest non matched sequence
     '''
    # we need the matches sorted by its start and merged when two overlap
    # we collect all match starts and ends in tuples
    locations = []
    for alignment in alignments:
        for match in alignment['matches']:
            for match_part in match['match_parts']:
                locations.append((match_part['query_start'],
                                  match_part['query_end']))
    # If no hsps, we return the complete sequence location
    if not locations:
        return (None, None)

    # Once we collect all the start and ends in a list of tuples, we are going
    # to merge overlapings
    locations = _merge_overlaping_locations(locations)

    # Now, we search for the longest non overlaping zone
    longest_location  = None
    longest_dist      = 0
    for i in range(len(locations)):
        if i == 0:
            start = 0
            end   = locations[i][0]
        else:
            start = locations[i-1][1]
            end   = locations[i][0]
        # Once we now the star and end of each non matching section we search
        # for the longest
        dist = end - start
        if dist > longest_dist:
            longest_dist = dist
            longest_location = (start, end)
    else:
        start  = locations[-1][1]
        end    = len(query)
        dist = end - start
        if dist > longest_dist:
            longest_location = (start, end)
    return longest_location

def create_masker_repeats_by_repeatmasker(species='eudicotyledons'):
    '''it creates a function that mask repeats giving only the seq, not the
    parameters'''

    def mask_repeats_by_repeatmasker(sequence):
        '''It returns a sequence with the repetitive and low complex elements
         masked'''
        if sequence is None:
            return None

        seq_fhand = run_repeatmasker_for_sequence(sequence, species)
        #pylint: disable-msg=W0612
        name, definition, masked_seq = get_content_from_fasta(seq_fhand)
        seq_class = sequence.seq.__class__
        return sequence.copy(seq=seq_class(masked_seq))
    return mask_repeats_by_repeatmasker

################################################################################
###  cleaner steps and specifications
#####################################

remove_vectors = {'function':create_vector_striper_by_alignment,
                  'arguments':{'vectors':None, 'aligner':'blast'},
                  'type': 'cleaner',
                  'name': 'remove_vectors',
                  'comment': 'Remove vector using vector db'}

remove_adaptors = {'function':create_vector_striper_by_alignment,
                  'arguments':{'vectors':None, 'aligner':'exonerate'},
                  'type': 'cleaner',
                  'name': 'remove_adaptors',
                  'comment': 'Remove our adaptors'}

strip_quality_lucy = {'function': create_striper_by_quality_lucy,
                      'arguments':{},
                      'type':'cleaner',
                      'name':'strip_lucy',
                      'comment':'Strip low quality with lucy'}

strip_quality_lucy2 = {'function': create_striper_by_quality_lucy2,
                      'arguments':{},
                      'type':'bulk_processor',
                      'name':'strip_lucy',
                      'comment':'Strip low quality with lucy'}

strip_quality_by_n = {'function': create_striper_by_quality_trimpoly,
                          'arguments': {},
                          'type':'cleaner',
                          'name':'strip_trimpoly',
                          'comment':'Strip low quality with trimpoly'}

mask_low_complexity = {'function': create_masker_for_low_complexity,
                       'arguments': {},
                       'type':'cleaner',
                       'name':'mask_low_complex',
                       'comment':'Mask low complexity regions'}

mask_repeats = {'function':create_masker_repeats_by_repeatmasker ,
                'arguments':{'species':'eudicotyledons'},
                'type': 'cleaner',
                'name': 'mask_repeats',
                'comment':'Mask repeats with repeatmasker'}

filter_short_seqs = {'function': create_length_filter,
                     'arguments':{'length':100, 'count_masked': False},
                     'type':'filter' ,
                     'name':'remove_short',
                     'comment': 'Remove seq shorter than 100 nt'}



PIPELINES = {'sanger_with_quality_clean'  : [remove_vectors, strip_quality_lucy2,
                                       mask_low_complexity, filter_short_seqs ],

            'sanger_without_quality_clean': [remove_vectors, strip_quality_by_n,
                                       mask_low_complexity, filter_short_seqs ],

            'repeatmasker'                :[mask_repeats, filter_short_seqs]
            }
def configure_pipeline(pipeline, configuration):
    '''It chooses the proper pipeline and configures it depending on the
    sequence kind and configuration parameters '''

    seq_pipeline  = PIPELINES[pipeline]

    # set the configuration in the pipeline
    for step in seq_pipeline:
        step_name = step['name']
        if step_name in configuration:
            for key, value in configuration[step_name].items():
                step['arguments'][key] = value

    # Here I check that none of the arguments have a none value
    for step in seq_pipeline:
        for key, value in step['arguments'].items():
            if value is None:
                msg = 'Parameter %s in step %s from pipeline %s must be set' % \
                            (key, step['name'], pipeline)
                raise RuntimeError(msg)
    return seq_pipeline


def pipeline_runner(pipeline, configuration, io_fhands, work_dir=None,
                    checkpoint_every_step=False):
    '''It runs all the analisis in the analisis_step especification dictionary

    This specidications depend on the type of the sequences'''

    # Here we extract our input/output files
    in_fhand_seqs  = io_fhands['in_seq']
    in_fhand_qual  = io_fhands['in_qual']
    out_fhand_seq  = io_fhands['out_seq']
    out_fhand_qual = io_fhands['out_qual']

    # We configure the pipeline depending on the sequences type and
    # configuratiom parameters
    pipeline_steps = configure_pipeline(pipeline, configuration)

    # Here starts the analisis
    seq_iter = seqs_in_file(in_fhand_seqs, in_fhand_qual)

    # List of temporary files created by the bulk processors.
    # we need to keep them until the analysis is done because some seq_iters
    # may depend on them
    temp_bulk_files = []

    for analisis_step in pipeline_steps:
        function_factory  = analisis_step['function']
        type_     = analisis_step['type']
        step_name = analisis_step['name']
        if analisis_step['arguments']:
            arguments = analisis_step['arguments']
        else:
            arguments = None

        msg = "Performing: %s" % analisis_step['comment']
        logging.info(msg)
        #print (msg)

        # Crete function adding parameters if they need them
        if arguments is None:
            cleaner_function = function_factory()
        else:
            # pylint:disable-msg=W0142
            cleaner_function = function_factory(**arguments)

        if type_ == 'cleaner':
            filtered_seqs = imap(cleaner_function, seq_iter)
        elif type_ == 'filter':
            filtered_seqs   = ifilter(cleaner_function, seq_iter)
        elif type_ == 'bulk_processor':
            filtered_seqs, fhand_outs = cleaner_function(seq_iter)
            temp_bulk_files.append(fhand_outs)

        # Now we need to create again the iterator. And this is going to be
        # useful to actually run the previous filter. Until now everything is
        # an iterator mega structure and nothing is executed
        if checkpoint_every_step:
            seq_step_name  = get_safe_fname(work_dir, step_name, 'seq.fasta')
            fhand_seq      = open(seq_step_name, 'w')
            if in_fhand_qual is not None:
                qual_step_name = get_safe_fname(work_dir, step_name,
                                                'qual.fasta')
                fhand_qual = open(qual_step_name, 'w')
            else:
                fhand_qual = None

            seq_iter = checkpoint(filtered_seqs, fhand_seq, fhand_qual)
            fhand_seq.close()
            if in_fhand_qual is not None:
                fhand_qual.close()
        else:
            seq_iter = filtered_seqs

        #more logging
        msg = "Finished: %s" % analisis_step['comment']
        logging.info(msg)

#        # I need an iterator to use with for the stats, so y generate another
#        # one with tee
#        seq_iter, seq_iter_stats = tee(seq_iter, 2)
#
#        #stats
#        run_stats(seq_iter_stats, analisis_step['statistics'], work_dir,
#                  step_name )

    else:
        checkpoint(seq_iter, out_fhand_seq, out_fhand_qual, 'false')

    # Cleaning saved temp files
    temp_bulk_files = None
    logging.info('Done!')

def run_stats(sequences, statistics, work_dir, step_name):
    '''It runs the statistics for the given sequences (iterator). Statistics is
    a list of stats to run'''
    if statistics is None:
        return
    msg = 'Running stats for %s' % step_name
    logging.info(msg)

    for statistic_func in statistics:
        sequences, new_sequences = tee(sequences, 2)
        stat_name = step_name + '.' + statistic_func.__name__

        distrib_fname = get_safe_fname(work_dir, stat_name, 'stat')
        plot_fname = get_safe_fname(work_dir, stat_name, 'png')

        statargs = {'sequences'     : sequences,
                    'distrib_fhand': open(distrib_fname, 'w'),
                    'plot_fhand'   : open(plot_fname, 'w')}
        # pylint:disable-msg=W0142
        statistic_func(**statargs)
        sequences = new_sequences

def checkpoint(seqs, out_fhand_seq, out_fhand_qual, return_iter=True):
    ''' This function is used to consume the previous iterators writing the
    output files. It generates another seq iterator that can be used for
    statistics or for the nerxt step.
    '''

    first      = True
    write_qual = None
    for seq in seqs:
        if seq is None:
            continue
        name     = seq.name
        sequence = seq.seq
        # check if we have quality and generate the fhand
        if first:
            first = False
            if seq.qual is None:
                write_qual = False
            else:
                write_qual = True
        # copy the seq to the fhand
        out_fhand_seq.write(fasta_str(sequence, name))
        if write_qual:
            quality  = ' '.join([str(qual) for qual in seq.qual])
            out_fhand_qual.write(fasta_str(quality, name))

    if return_iter:
        # files neew to be flushed
        out_fhand_seq.flush()
        if write_qual:
            out_fhand_qual.flush()
        # In order to be able t reas the generated file, It needs to be opened
        # for reading.
        seq_fname = out_fhand_seq.name
        fhand_seq_new = open(seq_fname, 'rt')

        if write_qual:
            qual_fname = out_fhand_qual.name
            fhand_qual_new = open(qual_fname, 'rt')
        else:
            fhand_qual_new = None

        return seqs_in_file(fhand_seq_new, fhand_qual_new)


