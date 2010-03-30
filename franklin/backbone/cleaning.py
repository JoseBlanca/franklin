'''
Created on 15/03/2010

@author: peio
'''
import os, itertools, logging
from franklin.backbone.analysis import Analyzer, scrape_info_from_fname
from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.backbone.specifications import BACKBONE_DIRECTORIES
from franklin.statistics import (seq_distrib, general_seq_statistics,
                                 seq_distrib_diff)
from franklin.seq.readers import guess_seq_file_format
from franklin.utils.seqio_utils import seqs_in_file

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

    def create_cleaning_configuration(self, platform, library=None):
        'It returns the pipeline configuration looking at the project settings'
        settings = self._project_settings['Cleaning']
        configuration = {}

        if 'vector_database' in settings:
            configuration['remove_vectors'] = {}
            configuration['remove_vectors']['vectors'] = \
                                                     settings['vector_database']

        # adaptors settings
        adaptors_file = None
        if platform == 'sanger' and 'adaptors_file_sanger' in settings:
            adaptors_file = settings['adaptors_file_sanger']
        elif platform == '454' and 'adaptors_file_454' in settings:
            adaptors_file = settings['adaptors_file_454']
        elif platform == 'illumina' and 'adaptors_file_illumina' in settings:
            adaptors_file = settings['adaptors_file_illumina']
        configuration['remove_adaptors'] = {}
        configuration['remove_adaptors']['vectors'] = adaptors_file


        # Words settings
        words = None
        if platform == 'sanger' and 'words_to_remove_sanger' in settings:
            words = settings['words_to_remove_sanger']
        elif platform == '454' and 'words_to_remove_454' in settings:
            words = settings['words_to_remove_454']
        elif platform == 'illumina' and 'words_to_remove_illumina' in settings:
            words = settings['words_to_remove_illumina']

        configuration['word_remover'] = {}
        configuration['word_remover']['words'] = words

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
        return configuration

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        self._log({'analysis_started':True})
        input_fpaths = self._get_input_fpaths()['reads']
        output_dir   =  self._create_output_dirs()['reads']
        for input_fpath in input_fpaths:
            fname = os.path.split(input_fpath)[-1]
            output_fpath = os.path.join(output_dir, fname)
            if os.path.exists(output_fpath):
                continue
            file_info = scrape_info_from_fname(input_fpath)
            input_fhand  = open(input_fpath)
            output_fhand  = open(output_fpath, 'w')
            pipeline = self._guess_cleaning_pipepile(file_info)
            iofhands = {'in_seq':input_fhand,
                        'outputs':{'sequence': output_fhand }}
            configuration = self.create_cleaning_configuration(
                                                       platform=file_info['pl'],
                                                       library=file_info['lb'])
            seq_pipeline_runner(pipeline, configuration, iofhands,
                                file_info['format'])
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

        clean_fpaths    = self._get_input_fpaths()['cleaned_reads']
        original_fpaths = self._get_input_fpaths()['original_reads']
        reads = {'cleaned': clean_fpaths, 'original': original_fpaths}

        # first per file stats
        for seq_type, fpaths in reads.items():
            for fpath in fpaths:
                stats_dir = self._get_stats_dir(seq_type)
                basename  = self._get_basename(fpath)
                seqs = seqs_in_file(open(fpath))
                file_format = guess_seq_file_format(open(fpath))
                analyses = ['seq_length_distrib', 'qual_distrib']
                if file_format == 'fasta':
                    analyses =  ['seq_length_distrib']

                self._do_seq_stats(seqs, basename, stats_dir, analyses)

        # now the difference between the cleaned and the original per file
        for clean_fpath in reads['cleaned']:
            original_fpath = os.path.join(self._get_project_path(),
                                         BACKBONE_DIRECTORIES['original_reads'],
                                         os.path.basename(clean_fpath))
            if os.path.exists(clean_fpath) and os.path.exists(original_fpath):
                stats_dir = self._get_stats_dir('cleaned')
                basename  = self._get_basename(clean_fpath)

                clean_seqs    = seqs_in_file(open(clean_fpath))
                original_seqs = seqs_in_file(open(original_fpath))
                file_format = guess_seq_file_format(open(clean_fpath))
                analyses = ['seq_length_distrib', 'qual_distrib']
                if file_format == 'fasta':
                    analyses =  ['seq_length_distrib']
                self._do_diff_seq_stats(clean_seqs, original_seqs, basename,
                                        stats_dir, analyses)

        # stats per seq file. All files together
        if clean_fpaths:
            clean_seqs = self._seqs_in_files(clean_fpaths)
            clean_seqs, clean_seqs2 = itertools.tee(clean_seqs, 2)
        else:
            clean_seqs, clean_seqs2 = None, None
        if original_fpaths:
            original_seqs = self._seqs_in_files(original_fpaths)
            original_seqs, original_seqs2 = itertools.tee(original_seqs, 2)
        else:
            original_seqs, original_seqs2 = None, None
        all_reads = {'cleaned': clean_seqs, 'original': original_seqs}
        basename  = 'global'
        for seqtype, seqs in all_reads.items():
            if seqs is None:
                continue
            stats_dir = self._get_stats_dir(seqtype)
            analyses = ['seq_length_distrib', 'qual_distrib']
            self._do_seq_stats(seqs, basename, stats_dir, analyses)

        # global stats. Diff fistributions
        stats_dir = self._get_stats_dir('cleaned')
        analyses  = ['seq_length_distrib', 'qual_distrib']
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
        plot_fpath    = os.path.join(stats_dir, basename + '.png')
        distrib_fpath = os.path.join(stats_dir, basename + '.dat')
        if os.path.exists(distrib_fpath) and os.path.exists(plot_fpath):
            return
        try:
            plot_fhand    = open(plot_fpath, 'w')
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
        plot_fpath    = os.path.join(stats_dir, basename + '.png')
        distrib_fpath = os.path.join(stats_dir, basename + '.dat')

        if os.path.exists(distrib_fpath) and os.path.exists(plot_fpath):
            return

        seq_distrib_diff(clean_seqs, original_seqs, analysis,
                         distrib_fhand=open(plot_fpath, 'w'),
                         plot_fhand=open(distrib_fpath, 'w'))

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
    'clean_read_stats':
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

