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

import os, itertools, logging
from tempfile import NamedTemporaryFile
from franklin.backbone.analysis import Analyzer, scrape_info_from_fname
from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.backbone.specifications import BACKBONE_DIRECTORIES
from franklin.statistics import (seq_distrib, general_seq_statistics,
                                 seq_distrib_diff)
from franklin.seq.readers import guess_seq_file_format
from franklin.utils.seqio_utils import seqs_in_file
from franklin.seq.seq_cleaner import MIN_LONG_ADAPTOR_LENGTH
from franklin.seq.writers import SequenceWriter

class CleanReadsAnalyzer(Analyzer):
    'It does a cleaning reads analysis'

    @staticmethod
    def _guess_cleaning_pipepile(file_info):
        'It returns the pipeline suited to clean the given file'
        if file_info['format'] == 'fasta':
            return 'sanger_without_qual'
        elif file_info['pl'] == 'illumina':
            return 'solexa'
        elif file_info['pl'] == '454':
            return 'sanger_with_qual'
        elif file_info['pl'] == 'sanger':
            return 'sanger_with_qual'
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

        if 'vector_database' in settings:
            configuration['remove_vectors'] = {}
            configuration['remove_vectors']['vectors'] = \
                                                     settings['vector_database']

        # adaptors settings
        adaptors_fpath = None
        adap_param = 'adaptors_file_%s' % platform
        if adap_param in settings:
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
        if word_param in settings:
            words.extend(settings[word_param])
        configuration['remove_short_adaptors'] = {}
        configuration['remove_short_adaptors']['words'] = words

        #edge_remover
        left, right = None, None
        if 'edge_removal' in settings:
            for (pl_side, length) in settings['edge_removal'].items():
                er_platform, side = pl_side.split('_')
                if er_platform == platform:
                    if side == 'left':
                        left = length
                    if side == 'right':
                        right = length
        configuration['edge_removal'] = {}
        configuration['edge_removal']['left_length'] = left
        configuration['edge_removal']['right_length'] = right

        # lucy settings.
        lucy_settings = settings['lucy_settings']
        if os.path.exists(lucy_settings):
            lucy_settings_dir = os.path.dirname(lucy_settings)
            lucy_settinhs_fhand = open(lucy_settings)
            lucy_libraries = eval(lucy_settinhs_fhand.read())
            lucy_settinhs_fhand.close()
            if library in lucy_libraries:
                vector = os.path.join(lucy_settings_dir,
                                      lucy_libraries[library]['vector_file'])
                splice = os.path.join(lucy_settings_dir,
                                      lucy_libraries[library]['splice_file'])
                configuration['strip_lucy'] = {}
                configuration['strip_lucy']['vector'] = [vector, splice]

        # min length settings
        min_seq_settings = settings['min_seq_length']
        if platform in min_seq_settings:
            min_length = min_seq_settings[platform]

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
            iofhands = {'in_seq':input_fhand,
                        'outputs':{'sequence': output_fhand }}
            configuration = self.create_cleaning_configuration(
                                                       platform=file_info['pl'],
                                                       library=file_info['lb'])
            seq_pipeline_runner(pipeline, configuration, iofhands,
                                file_info['format'], processes=self.threads)
            input_fhand.close()
            output_fhand.close()
        self._log({'analysis_finished':True})
        return

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
        reads = {'cleaned': clean_paths, 'original': original_paths}

        if ('reads_stats' in self._project_settings and
            'sampling_size' in self._project_settings['reads_stats' ]):
            sample_size = self._project_settings['reads_stats' ]['sampling_size']
        else:
            sample_size = None

        # first per file stats
        for seq_type, paths in reads.items():
            for path in paths:
                stats_dir = self._get_stats_dir(seq_type)
                basename = path.basename
                fpath = path.last_version
                seqs = seqs_in_file(open(fpath), sample_size=sample_size)
                file_format = guess_seq_file_format(open(fpath))
                analyses = ['seq_length_distrib', 'qual_distrib']
                if file_format == 'fasta':
                    analyses = ['seq_length_distrib']

                self._do_seq_stats(seqs, basename, stats_dir, analyses)

        # now the difference between the cleaned and the original per file
        for clean_path in reads['cleaned']:
            original_fpath = os.path.join(self._get_project_path(),
                                         BACKBONE_DIRECTORIES['original_reads'],
                                         clean_path.basename)
            clean_fpath = clean_path.last_version
            if os.path.exists(clean_fpath) and os.path.exists(original_fpath):
                stats_dir = self._get_stats_dir('cleaned')

                clean_seqs = seqs_in_file(open(clean_fpath),
                                          sample_size=sample_size)
                original_seqs = seqs_in_file(open(original_fpath),
                                             sample_size=sample_size)
                file_format = guess_seq_file_format(open(clean_fpath))
                analyses = ['seq_length_distrib', 'qual_distrib']
                if file_format == 'fasta':
                    analyses = ['seq_length_distrib']
                self._do_diff_seq_stats(clean_seqs, original_seqs,
                                        clean_path.basename,
                                        stats_dir, analyses)

        # stats per seq file. All files together
        if clean_paths:
            clean_fpaths = [path.last_version for path in clean_paths]
            clean_seqs = self._seqs_in_files(clean_fpaths)
            clean_seqs, clean_seqs2 = itertools.tee(clean_seqs, 2)
        else:
            clean_seqs, clean_seqs2 = None, None
        if original_paths:
            original_fpaths = [path.last_version for path in original_paths]
            original_seqs = self._seqs_in_files(original_fpaths)
            original_seqs, original_seqs2 = itertools.tee(original_seqs, 2)
        else:
            original_seqs, original_seqs2 = None, None
        all_reads = {'cleaned': clean_seqs, 'original': original_seqs}
        basename = 'global'
        for seqtype, seqs in all_reads.items():
            if seqs is None:
                continue
            stats_dir = self._get_stats_dir(seqtype)
            analyses = ['seq_length_distrib', 'qual_distrib']
            self._do_seq_stats(seqs, basename, stats_dir, analyses)

        # global stats. Diff fistributions
        stats_dir = self._get_stats_dir('cleaned')
        analyses = ['seq_length_distrib', 'qual_distrib']
        if clean_seqs2 is not None and original_seqs2 is not None:
            self._do_diff_seq_stats(clean_seqs2, original_seqs2, basename,
                                    stats_dir, analyses)
        self._log({'analysis_finished':True})

    def _get_stats_dir(self, seqtype):
        'It gets the stats dir for each seqtype'
        return os.path.join(self._get_project_path(),
                             BACKBONE_DIRECTORIES['%s_reads_stats' % seqtype])
    @staticmethod
    def _seqs_in_files(fpaths):
        'It yields seqrecored from a list of files'
        for fpath in fpaths:
            for seq in seqs_in_file(open(fpath)):
                yield seq

    def _do_seq_stats(self, seqs, basename, stats_dir, analyses):
        'It performs all kind of stats for a fpath'
        # Some distributions
        for analysis in analyses:
            analysis_basename = '%s.%s' % (basename, analysis)
            seqs, seqs2 = itertools.tee(seqs, 2)
            self._do_distrib(seqs2, analysis, analysis_basename, stats_dir)

        # now general stats
        analysis_basename = '%s.general_stats' % basename
        self._do_general_stats(seqs, analysis_basename, stats_dir)

    def _do_general_stats(self, seqs, basename, stats_dir):
        'It performs the  general stats analysis'
        general_fpath = os.path.join(self._get_project_path(), stats_dir,
                                     basename + '.dat')
        if os.path.exists(general_fpath):
            return

        general_fhand = open(general_fpath, 'w')
        stats = general_seq_statistics(seqs)
        for key, value in stats.items():
            if value is not None:
                to_print = '%-19s : %d\n' % (key, value)
                general_fhand.write(to_print)

    @staticmethod
    def _do_distrib(seqs, analysis, basename, stats_dir):
        'I actually do the distrib'
        plot_fpath = os.path.join(stats_dir, basename + '.png')
        distrib_fpath = os.path.join(stats_dir, basename + '.dat')
        if os.path.exists(distrib_fpath) and os.path.exists(plot_fpath):
            return
        try:
            plot_fhand = open(plot_fpath, 'w')
            distrib_fhand = open(distrib_fpath, 'w')

            seq_distrib(sequences=seqs, kind=analysis,
                        distrib_fhand=distrib_fhand, plot_fhand=plot_fhand,
                        low_memory=True)
        except TypeError:
            plot_fhand.close()
            distrib_fhand.close()
            os.remove(plot_fpath)
            os.remove(distrib_fpath)
            raise

    def _do_diff_seq_stats(self, seqs1, seqs2, basename, stats_dir, analyses):
        'It performs the differential distribution'

        for analysis in analyses:
            seqs1, seqs1_use = itertools.tee(seqs1, 2)
            seqs2, seqs2_use = itertools.tee(seqs2, 2)
            analysis_basename = '%s.diff_%s' % (basename, analysis)
            self._do_diff_distrib(seqs1_use, seqs2_use, analysis,
                                  analysis_basename, stats_dir)

    @staticmethod
    def _do_diff_distrib(clean_seqs, original_seqs, analysis, basename,
                         stats_dir):
        'It performs the differential distribution'
        plot_fpath = os.path.join(stats_dir, basename + '.png')
        distrib_fpath = os.path.join(stats_dir, basename + '.dat')

        if os.path.exists(distrib_fpath) and os.path.exists(plot_fpath):
            return

        seq_distrib_diff(clean_seqs, original_seqs, analysis,
                         distrib_fhand=open(distrib_fpath, 'w'),
                         plot_fhand=open(plot_fpath, 'w'))

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
         'outputs':{'original_reads':{'directory':'original_reads_stats'},
                    'cleaned_reads'   :{'directory':'cleaned_reads_stats'}},
         'analyzer': ReadsStatsAnalyzer,
        },
     }
