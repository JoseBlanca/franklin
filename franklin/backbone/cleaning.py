'''
This module is part of ngs_backbone. This module performs analyses related to
sequence cleaning.

Created on 15/03/2010

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

import os, logging, array
from tempfile import NamedTemporaryFile

from franklin.backbone.analysis import Analyzer, scrape_info_from_fname
from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.backbone.specifications import (BACKBONE_DIRECTORIES,
                                              BACKBONE_BASENAMES,
                                              PLOT_FILE_FORMAT)
from franklin.utils.seqio_utils import seqs_in_file
from franklin.seq.seq_cleaner import MIN_LONG_ADAPTOR_LENGTH
from franklin.seq.writers import SequenceWriter
from franklin.statistics import (create_distribution, write_distribution,
                                 draw_histogram, CachedArray, draw_boxplot)

class CleanReadsAnalyzer(Analyzer):
    'It does a cleaning reads analysis'

    @staticmethod
    def _guess_cleaning_pipepile(file_info):
        'It returns the pipeline suited to clean the given file'
        if 'pl' not in file_info:
            raise RuntimeError('platform not in the file name: %s' %
                               file_info['fpath'])
        platform = file_info['pl']
        if file_info['format'] == 'fasta':
            return 'sanger_without_qual'
        elif platform == 'illumina':
            return 'solexa'
        elif platform == '454':
            return 'sanger_with_qual'
        elif platform == 'sanger':
            return 'sanger_with_qual'
        elif platform == 'solid':
            return 'solid'
        else:
            raise ValueError('Unable to guess the cleaning pipeline: %s' %
                             str(file_info))

    @staticmethod
    def _purge_short_adaptors(adaptors_fpath):
        '''It returns a file without short adaptors and the list of adaptors
        removed'''
        long_adap_fhand = NamedTemporaryFile(delete=False, mode='a')
        seq_writer = SequenceWriter(long_adap_fhand, file_format='fasta')
        short_adaptors = []
        for seq in seqs_in_file(open(adaptors_fpath)):
            if len(seq) < MIN_LONG_ADAPTOR_LENGTH:
                short_adaptors.append(str(seq.seq))
            else:
                seq_writer.write(seq)
        if not short_adaptors:
            #not worth the new file, we remove it
            long_adap_fhand.close()
            os.remove(long_adap_fhand.name)
            return adaptors_fpath, short_adaptors
        elif not seq_writer.num_features:   #no long adaptors remaining
            os.remove(long_adap_fhand.name)
            return None, short_adaptors     #we do not return the adaptor file
        else:
            return long_adap_fhand.name, short_adaptors

    def create_cleaning_configuration(self, platform, library=None):
        'It returns the pipeline configuration looking at the project settings'
        settings = self._project_settings['Cleaning']
        configuration = {}

        configuration['remove_vectors'] = {}
        configuration['remove_vectors']['vectors'] = settings['vector_database']

        # adaptors settings
        adap_param = 'adaptors_file_%s' % platform
        adaptors_fpath = settings[adap_param]

        #we have to remove the short adaptors from the file and treat them
        #as short adaptors
        if adaptors_fpath:
            adaptors_fpath, words = self._purge_short_adaptors(adaptors_fpath)
        else:
            words = []
        configuration['remove_adaptors'] = {}
        configuration['remove_adaptors']['vectors'] = adaptors_fpath

        # Words settings
        word_param = 'short_adaptors_%s' % platform
        if settings[word_param] is not None:
            words.extend(settings[word_param])
        configuration['remove_short_adaptors'] = {}
        configuration['remove_short_adaptors']['words'] = words

        #edge_remover
        left =  settings['edge_removal']['%s_left' % platform]
        right = settings['edge_removal']['%s_right' % platform]
        configuration['edge_removal'] = {}
        configuration['edge_removal']['left_length'] = left
        configuration['edge_removal']['right_length'] = right

        # lucy settings.
        lucy_settings = settings['lucy']
        lucy_vector_settings = lucy_settings['vector_settings']
        if lucy_vector_settings is not None:
            lucy_settings_dir = os.path.dirname(lucy_vector_settings)
            lucy_settings_fhand = open(lucy_vector_settings)
            lucy_libraries = eval(lucy_settings_fhand.read())
            lucy_settings_fhand.close()
            if library in lucy_libraries:
                vector = os.path.join(lucy_settings_dir,
                                      lucy_libraries[library]['vector_file'])
                splice = os.path.join(lucy_settings_dir,
                                      lucy_libraries[library]['splice_file'])
                vector_settings = [vector, splice]
            else:
                vector_settings = None
        else:
            vector_settings = None

        parameters = {'bracket':lucy_settings['bracket'],
                      'window':lucy_settings['window'],
                      'error':lucy_settings['error']}
        if vector_settings:
            parameters['vector'] = vector_settings

        configuration['strip_lucy'] = {}
        configuration['strip_lucy']['parameters'] = parameters

        n_trim = settings['strip_n_percent']
        configuration['strip_trimpoly'] = {}
        configuration['strip_trimpoly']['ntrim_above_percent'] = n_trim

        # min length settings
        min_length = settings['min_seq_length'][platform]
        configuration['remove_short'] = {}
        configuration['remove_short']['length'] = min_length

        return configuration

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        logger = logging.getLogger("franklin")
        self._log({'analysis_started':True})
        input_paths = self._get_input_fpaths()['reads']
        output_dir = self._create_output_dirs()['reads']
        for input_path in input_paths:
            input_fpath = str(input_path)
            fname = os.path.split(input_fpath)[-1]
            output_fpath = os.path.join(output_dir, fname)
            if os.path.exists(output_fpath):
                msg = '%s already cleaned. Not cleaned again'  % output_fpath
                logger.info(msg)
                continue
            file_info = scrape_info_from_fname(input_path)
            input_fhand = open(input_fpath)
            output_fhand = open(output_fpath, 'w')
            pipeline = self._guess_cleaning_pipepile(file_info)
            infhands = {'in_seq':input_fhand}
            writer = SequenceWriter(output_fhand,
                                    file_format=file_info['format'])

            configuration = self.create_cleaning_configuration(
                                                       platform=file_info['pl'],
                                                       library=file_info['lb'])
            seq_pipeline_runner(pipeline, configuration, infhands,
                                file_info['format'], processes=self.threads,
                                writers={'seq':writer})
            input_fhand.close()
            output_fhand.close()
        self._log({'analysis_finished':True})
        return

BASENAME_FOR_ALL_TOGETHER = 'all'

PLOT_LABELS = {
        'seq_length':{'title':'Sequence length distribution',
                              'xlabel':'Sequence length',
                              'ylabel': 'Number of sequences'},
        'seq_qual'      :{'title':'Sequence quality distribution',
                              'xlabel':'quality',
                              'ylabel':'Number of bp'},
        'qual_boxplot'  :{'title':'Qualities along the sequence',
                              'xlabel':'seqence position',
                              'ylabel':'quality distribution'},
        }

class ReadsStatsAnalyzer(Analyzer):
    '''It calculates stats for original and cleaned reads.
    It calculates the distribution between both type off reads'''

    def run(self):
        'It runs the analysis'
        self._log({'analysis_started':True})
        self._create_output_dirs()['original_reads']
        self._create_output_dirs()['cleaned_reads']

        clean_paths = self._get_input_fpaths()['cleaned_reads']
        original_paths = self._get_input_fpaths()['original_reads']
        reads = {'cleaned': clean_paths, 'raw': original_paths}

        #we need the raw and cleaned files paired
        paired_paths = self._pair_raw_cleaned_read_files(reads)

        #now we can do the statistics
        for pair in paired_paths.values():
            self._do_seq_distrib_for_pair(pair)

        self._log({'analysis_finished':True})

    @staticmethod
    def _do_diff_distrib_for_numbers(numbers, plot_fhand, distrib_fhand,
                                     dist_type):
        'It creates the diff distribution for the given numbers'

        max_ = max(numbers[0].max, numbers[1].max)
        min_ = min(numbers[0].min, numbers[1].min)
        range_ = min_, max_
        # to get the difference we need both distribs
        distrib1 = create_distribution(numbers[0], range_=range_)
        distrib2 = create_distribution(numbers[1], range_=range_)

        # now a subtract distrib1 from distrib2
        diff_distrib   = []
        diff_bin_edges = distrib1['bin_edges']
        for i in range(len(distrib1['distrib'])):
            diff = distrib1['distrib'][i] - distrib2['distrib'][i]
            diff_distrib.append(diff)

        #now we write the result
        write_distribution(distrib_fhand, diff_distrib, diff_bin_edges)

        labels = PLOT_LABELS[dist_type]
        draw_histogram(diff_distrib, diff_bin_edges,
                       title=labels['title'], xlabel=labels['xlabel'],
                       ylabel=labels['ylabel'],
                       fhand=plot_fhand)

    def _do_seq_distrib_for_pair(self, pair):
        'It does the distribution for a pair of cleaned and raw seqs'

        get_stats_dir = lambda seq_type: os.path.join(self._get_project_path(),
                              BACKBONE_DIRECTORIES['%s_reads_stats' % seq_type])

        #the statistics for both clean and raw sequences
        lengths = {}
        quals = {}
        for seq_type in ('raw', 'cleaned'):
            stats_dir = get_stats_dir(seq_type)

            stats_fpath = os.path.join(stats_dir,
                                       BACKBONE_BASENAMES['statistics_file'])
            stats_fhand = open(stats_fpath, 'a')

            if seq_type in pair:
                fpath = pair[seq_type].last_version
                basename = pair[seq_type].basename

                #the names for the output files
                out_fpath = os.path.join(stats_dir, basename + '.length')
                plot_fpath = out_fpath + '.' + PLOT_FILE_FORMAT
                distrib_fpath = out_fpath + '.dat'

                if os.path.exists(plot_fpath):
                    continue

                lengths_, quals_ = self._get_lengths_quals_from_file(fpath)
                lengths[seq_type] = lengths_
                quals[seq_type] = quals_

                #the distributions for the lengths
                create_distribution(lengths_, PLOT_LABELS['seq_length'],
                                    distrib_fhand=open(distrib_fpath, 'w'),
                                    plot_fhand=open(plot_fpath, 'w'),
                                    range_= (lengths_.min, lengths_.max))

                #the distributions for the quals
                out_fpath = os.path.join(stats_dir, basename + '.qual')
                plot_fpath = out_fpath + '.' + PLOT_FILE_FORMAT
                distrib_fpath = out_fpath + '.dat'
                if quals_:
                    create_distribution(quals_, PLOT_LABELS['seq_qual'],
                                        distrib_fhand=open(distrib_fpath, 'w'),
                                        plot_fhand=open(plot_fpath, 'w'),
                                        range_=(quals_.min, quals_.max))

                #the statistics for the statistics file
                self._write_statistics(stats_fhand, fpath, lengths_, quals_)

        #the statistics for the differences
        if 'raw' in pair and 'cleaned' in pair:
            fpath = pair['cleaned'].last_version
            basename = pair['cleaned'].basename

            #the names for the output files
            stats_dir = get_stats_dir('cleaned')
            out_fpath = os.path.join(stats_dir, basename + '.length.diff')
            plot_fpath = out_fpath + '.' + PLOT_FILE_FORMAT
            distrib_fpath = out_fpath + '.dat'

            if not os.path.exists(plot_fpath):
                #the distributions for the lengths
                lengths = lengths['raw'], lengths['cleaned']
                self._do_diff_distrib_for_numbers(lengths,
                                          plot_fhand= open(plot_fpath, 'w'),
                                    distrib_fhand= open(distrib_fpath, 'w'),
                                                  dist_type='seq_length')
                del lengths

                #the distributions for the quals
                out_fpath = os.path.join(stats_dir, basename + '.qual.diff')
                plot_fpath = out_fpath + '.' + PLOT_FILE_FORMAT
                distrib_fpath = out_fpath + '.dat'

                quals = quals['raw'], quals['cleaned']
                if quals[0]:
                    self._do_diff_distrib_for_numbers(quals,
                                          plot_fhand= open(plot_fpath, 'w'),
                                    distrib_fhand= open(distrib_fpath, 'w'),
                                                      dist_type='seq_qual')
                del quals

        for seq_type in ('raw', 'cleaned'):
            if seq_type in pair:
                fpath = pair[seq_type].last_version
                basename = pair[seq_type].basename

                #the names for the output files
                out_fpath = os.path.join(stats_dir, basename + '.qual.boxplot')
                plot_fpath = out_fpath + '.' + PLOT_FILE_FORMAT
                distrib_fpath = out_fpath + '.dat'

                if os.path.exists(plot_fpath):
                    continue

                quals_ = self._get_quals_by_length_from_file(fpath)

                if quals_ and quals_[0]:
                    #boxplot
                    draw_boxplot(quals_, fhand=open(plot_fpath, 'w'),
                                 title=PLOT_LABELS['qual_boxplot']['title'],
                                 xlabel=PLOT_LABELS['qual_boxplot']['xlabel'],
                                 ylabel=PLOT_LABELS['qual_boxplot']['ylabel'],
                                 stats_fhand=open(distrib_fpath, 'w'),
                                 max_plotted_boxes=30)
                del quals_

    @staticmethod
    def _get_quals_by_length_from_file(fpath):
        'It returns the qualities along the sequence as a list of lists'

        quals_by_position = []
        for seq in seqs_in_file(open(fpath)):
            quals = seq.qual
            if quals:
                for base_number, qual in enumerate(quals):
                    try:
                        quals_by_position[base_number]
                    except IndexError:
                        quals_by_position.append(array.array('B'))
                    quals_by_position[base_number].append(qual)
        return quals_by_position

    @staticmethod
    def _write_statistics(stats_fhand, seq_fpath, lengths, quals):
        'It writes some statistics'

        to_print = 'statistics for %s\n' % os.path.basename(seq_fpath)
        stats_fhand.write(to_print)
        stats_fhand.write('-' * (len(to_print) - 1) + '\n')

        stats_fhand.write('Num sequences: %i\n' % len(lengths))

        stats_fhand.write('Total sequence length: %i\n' % lengths.sum)

        stats_fhand.write('Sequence length minimum: %i\n' % lengths.min)

        stats_fhand.write('Sequence length maximum: %i\n' % lengths.max)

        stats_fhand.write('Sequence length average: %.2f\n' % lengths.average)

        stats_fhand.write('Sequence length variance: %.2f\n' % lengths.variance)
        if quals:
            stats_fhand.write('Sequence qualities minimum: %i\n' % quals.min)

            stats_fhand.write('Sequence qualities maximum: %i\n' % quals.max)

            stats_fhand.write('Sequence qualities average: %.2f\n' % \
                                                                  quals.average)
            stats_fhand.write('Sequence qualities variance: %.2f\n' % \
                                                                 quals.variance)
        stats_fhand.write('\n')

    @staticmethod
    def _get_lengths_quals_from_file(seq_fpath):
        'Given a sequence file it returns the lengths and quals'
        lengths = CachedArray('I')
        quals   = CachedArray('H')
        for seq in seqs_in_file(open(seq_fpath)):
            lengths.append(len(seq))
            qual = seq.qual
            if qual:
                quals.extend(qual)
        return lengths, quals

    @staticmethod
    def _pair_raw_cleaned_read_files(read_paths):
        'It returns the paths for cleaned and raw paired'
        pairs = {}
        for path in read_paths['raw']:
            basename = path.basename
            pairs[basename] = {'raw': path}
        for path in read_paths['cleaned']:
            basename = path.basename
            if basename not in pairs:
                pairs[basename] = {}
            pairs[basename]['cleaned'] = path
        return pairs

DEFINITIONS = {
    'clean_reads':
        {'inputs':{
            'reads':
                {'directory': 'original_reads',
                 'file_kinds': 'sequence_files'}
            },
         'outputs':{'reads':{'directory': 'cleaned_reads'}},
         'analyzer': CleanReadsAnalyzer,
        },
    'read_stats':
        {'inputs':{
            'original_reads':
                {'directory': 'original_reads',
                 'file_kinds': 'sequence_files'},
            'cleaned_reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'}
                },
         'outputs':{'original_reads':{'directory':'raw_reads_stats'},
                    'cleaned_reads'   :{'directory':'cleaned_reads_stats'}},
         'analyzer': ReadsStatsAnalyzer,
        },
     }
