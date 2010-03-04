'''
Created on 01/03/2010

@author: peio
'''


import os, time, shutil, itertools, tempfile
from configobj import ConfigObj
from tempfile import NamedTemporaryFile
from franklin.utils.cmd_utils import call
from franklin.utils.misc_utils import NamedTemporaryDir
from franklin.utils.seqio_utils import (seqio, cat, seqs_in_file,
                                        write_seqs_in_file)
from franklin.seq.readers import guess_seq_file_format
from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.mapping import map_reads
from franklin.sam import (bam2sam, add_header_and_tags_to_sam, merge_sam,
                          sam2bam, sort_bam_sam)
from franklin.snv.sam_pileup import snv_contexts_in_sam_pileup
from franklin.pipelines.pipelines import _pipeline_builder
from franklin.statistics import (seq_distrib, general_seq_statistics,
                                 seq_distrib_diff)


def _is_file_kind(fpath, extensions):
    'It returns True if the file has one of the given extensions'
    ext = os.path.splitext(fpath)[-1].strip('.')
    if ext in extensions:
        return True
    else:
        return False

def _is_sequence_file(fpath):
    'It returns true if the function is a sequence'
    return _is_file_kind(fpath,  ('fasta', 'fastq', 'sfastq', 'repr'))

def _is_sequence_or_qual_file(fpath):
    'It returns true if the function is a sequence or quality file'
    return _is_file_kind(fpath,  ('fasta', 'fastq', 'sfastq', 'repr', 'qual'))

def _select_file(kind, fpaths):
    'It returns the fpath that correpond to the given file kind'
    if not fpaths:
        raise RuntimeError('No fpaths given, so no selection possible')
    if kind == 'contigs':
        fpaths = filter(_is_sequence_file, fpaths)
        fpaths = filter(lambda x: BACKBONE_BASENAMES['contigs'] in x, fpaths)
    elif  kind == 'mapping_reference':
        fpaths = filter(_is_sequence_file, fpaths)
        fpaths = filter(lambda x: BACKBONE_BASENAMES['mapping_reference'] in x,
                        fpaths)
    elif kind == 'merged_bam':
        fpaths = filter(lambda x: BACKBONE_BASENAMES['merged_bam'] in x,
                        fpaths)
    if not fpaths:
        raise RuntimeError('No file found for %s' % kind)
    assert len(fpaths) == 1
    return fpaths[0]

class Analyzer(object):
    '''This class performs an analysis.

    This is a prototype that should be inherited for each analysis. This class
    looks for inputs files and runs the analysis.
    It checks if the the analysis is already done, looking to the output
    directory
    '''
    def __init__(self, project_settings, analysis_definition):
        'The init'
        self._project_settings = project_settings
        self._analysis_def = analysis_definition
        self._timestamped_dir = None
        self._old_tmpdir = tempfile.gettempdir()
        self._setup_tempdir()

    @staticmethod
    def _set_tmp(tmpdir):
        'It sets the tmpdir'
        if tmpdir is None:
            raise NotImplementedError
        else:
            if not os.path.exists(tmpdir):
                os.makedirs(tmpdir)
            tempfile.tempdir = tmpdir
            os.environ['TMPDIR'] = tmpdir

    def _setup_tempdir(self):
        'It prepares tempfile to use the given tempdir from the settings'

        self._old_tempdir = tempfile.tempdir
        try:
            tmpdir = self._project_settings['General_settings']['tmpdir']
        except KeyError:
            return
        self._set_tmp(tmpdir)

    def __del__(self):
        'Some clean up'
        self._set_tmp(self._old_tmpdir)

    def _get_project_name(self):
        'It returns the name of the project'
        return self._project_settings['General_settings']['project_name']

    def _get_project_path(self):
        'It return the project path'
        return self._project_settings['General_settings']['project_path']

    def _get_input_dirs(self):
        'It returns the directories for the inputs'
        return self._get_io_dirs('inputs')

    def _get_output_dirs(self, timestamped=False):
        'It returns the directories for the output'
        return self._get_io_dirs('outputs', timestamped)

    def _get_io_dirs(self, stream_kind, timestamped=False):
        'It returns the inputs or outputs directories'
        stream_defs = self._analysis_def[stream_kind]
        dirs = {}
        project_dir = self._project_settings['General_settings']['project_path']
        for kind, stream_def in stream_defs.items():
            output_dir = BACKBONE_DIRECTORIES[stream_def['directory']]
            if timestamped:
                if self._timestamped_dir:
                    time_bit = self._timestamped_dir
                else:
                    time_bit = time.strftime("%Y%m%d_%H%M", time.gmtime())
                    self._timestamped_dir = time_bit
            else:
                time_bit = ''
            if isinstance(output_dir, tuple):
                pre_time_bit, post_time_bit = output_dir
            else:
                pre_time_bit = output_dir
                post_time_bit = ''

            timed_dir = os.path.join(project_dir,
                                     pre_time_bit,
                                     time_bit,
                                     post_time_bit)
            dirs[kind] = timed_dir
        return dirs

    def _get_input_fpaths(self):
        'It returns a dict with the inputs files'
        inputs_def = self._analysis_def['inputs']
        input_files = {}
        input_dirs = self._get_input_dirs()
        for input_kind, input_def in inputs_def.items():
            backbone_dir = input_dirs[input_kind]
            fpaths = os.listdir(backbone_dir)

            #we want full paths
            fpaths = [os.path.join(backbone_dir, path) for path in fpaths]

            if 'file_kinds' in input_def:
                if input_def['file_kinds'] == 'sequence_files':
                    fpaths = filter(_is_sequence_file, fpaths)
                elif input_def['file_kinds'] == 'sequence_or_qual_files':
                    fpaths = filter(_is_sequence_or_qual_file, fpaths)
                elif input_def['file_kinds'] == 'bam':
                    fpaths = filter(lambda x: x.endswith('.bam'), fpaths)
                elif input_def['file_kinds'] == 'pileup':
                    fpaths = filter(lambda x: x.endswith('.pileup'), fpaths)
                elif input_def['file_kinds'] == 'frg_file':
                    fpaths = filter(lambda x: x.endswith('.frg'), fpaths)
            if 'file' in input_def:
                fpaths = _select_file(input_def['file'], fpaths)
            input_files[input_kind] = fpaths
        return input_files

    def get_output_fpaths(self):
        'It returns the output files for the analysis'
        raise NotImplementedError

    def run(self):
        'It runs the analysis'
        raise NotImplementedError

    def _get_timestamped_output_dirs(self):
        'It returns the path for the output dir in each timestamped result'
        analysis_dir = self._get_input_dirs()['analyses_dir']

        def _is_timestamped_dir(fname):
            'It returns true for the timestamped directories'
            if fname[0] == '2'and fname[-1].isdigit():
                return True
            else:
                return False
        dir_names = filter(_is_timestamped_dir, os.listdir(analysis_dir))
        #for each timestamped dir there is a /result directory, this should
        #be part of the dir_name
        result_dir_kind = self._analysis_def['outputs']['result']['directory']
        result_dir = BACKBONE_DIRECTORIES[result_dir_kind][1]
        dir_names = [os.path.join(dname, result_dir) for dname in dir_names]

        dirs = [os.path.join(analysis_dir, dname) for dname in dir_names]
        return sorted(dirs)


    def _create_output_dirs(self, timestamped=False):
        'It creates the output directories if they do not exists'
        output_dirs = self._get_output_dirs(timestamped)
        for kind, output_dir in output_dirs.items():
            if ('create' in self._analysis_def['outputs'][kind] and
                not self._analysis_def['outputs'][kind]['create']):
                continue
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        return output_dirs

class CleanReadsAnalyzer(Analyzer):
    'It does a cleaning reads analysis'

    def _guess_cleaning_pipepile(self, file_info):
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

        adaptors_file = None
        if platform == 'sanger' and 'adaptors_file_sanger' in settings:
            adaptors_file = settings['adaptors_file_sanger']
        elif platform == '454' and 'adaptors_file_454' in settings:
            adaptors_file = settings['adaptors_file_454']
        elif platform == 'illumina' and 'adaptors_file_illumina' in settings:
            adaptors_file = settings['adaptors_file_illumina']
        if adaptors_file:
            configuration['remove_adaptors'] = {}
            configuration['remove_adaptors']['vectors'] = adaptors_file

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
        input_fpaths = self._get_input_fpaths()['reads']
        output_dir   =  self._create_output_dirs()['reads']
        for input_fpath in input_fpaths:
            fname = os.path.split(input_fpath)[-1]
            output_fpath = os.path.join(output_dir, fname)
            if os.path.exists(output_fpath):
                continue
            file_info = _scrape_info_from_fname(input_fpath)
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
        return

def _scrape_info_from_fname(fpath):
    'It guess pipeline taking into account the platform and the file format'
    file_info = {}
    fhand = open(fpath)
    file_info['format'] = guess_seq_file_format(fhand)
    fhand.close()
    fname = os.path.splitext(os.path.split(fpath)[-1])[0]
    for item in fname.split('.'):
        key, value     = item.split('_', 1)
        file_info[key] = value
    return file_info


class ReadsStatsAnalyzer(Analyzer):
    '''It calculates stats for original and cleaned reads.
    It calculates the distribution between both type off reads'''

    def run(self):
        'It runs the analysis'
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
        clean_seqs    = self._seqs_in_files(clean_fpaths)
        original_seqs = self._seqs_in_files(original_fpaths)
        clean_seqs, clean_seqs2       = itertools.tee(clean_seqs, 2)
        original_seqs, original_seqs2 = itertools.tee(original_seqs, 2)

        all_reads = {'cleaned': clean_seqs, 'original': original_seqs}
        basename  = 'global'
        for seqtype, seqs in all_reads.items():
            stats_dir = self._get_stats_dir(seqtype)
            analyses = ['seq_length_distrib', 'qual_distrib']
            self._do_seq_stats(seqs, basename, stats_dir, analyses)

        # global stats. Diff fistributions
        stats_dir = self._get_stats_dir('cleaned')
        self._do_diff_seq_stats(clean_seqs2, original_seqs2, basename,
                                stats_dir)

    @staticmethod
    def _get_basename(fpath):
        'It gets the basename of the file and the stats'
        return os.path.splitext(os.path.basename(fpath))[0]

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

    def _do_distrib(self, seqs, analysis, basename, stats_dir):
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

    def _do_diff_distrib(self, clean_seqs, original_seqs, analysis, basename,
                         stats_dir):
        'It performs the differential distribution'
        plot_fpath    = os.path.join(stats_dir, basename + '.png')
        distrib_fpath = os.path.join(stats_dir, basename + '.dat')

        if os.path.exists(distrib_fpath) and os.path.exists(plot_fpath):
            return

        seq_distrib_diff(clean_seqs, original_seqs, analysis,
                         distrib_fhand=open(plot_fpath, 'w'),
                         plot_fhand=open(distrib_fpath, 'w'))


class PrepareWSGAssemblyAnalyzer(Analyzer):
    '''It collects the cleaned reads to use by wsg. Wsg only uses reads with quality,
    so be will dismiss these sequences'''

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''

        # TODO
        # Now we need to convert all the fastq in fasta and qual. If fasta
        #without qual, we need to give it a quality from settings. It can change
        # if people in celera get fastqToCA to work.
        tempdir = NamedTemporaryDir()
        tempdir_name = tempdir.name
        fasta_qual_files = []
        for seqfile in self._get_input_fpaths()['reads']:
            seqfhand = open(seqfile)
            basename = os.path.splitext(os.path.basename(seqfile))[0]
            fasta_fhand = open(os.path.join(tempdir_name,
                                            basename + '.fasta'), 'w')
            qual_fhand  = open(os.path.join(tempdir_name,
                                            basename + '.qual'), 'w')
            fasta_qual_files.append((fasta_fhand.name, qual_fhand.name))

            file_format = guess_seq_file_format(seqfhand)
            seqs = seqs_in_file(seqfhand, format=file_format)

            write_seqs_in_file(seqs,  seq_fhand=fasta_fhand,
                               qual_fhand=qual_fhand,
                               format='fasta')
            fasta_fhand.close()
            qual_fhand.close()
            seqfhand.close()

        # once we have the fasta and qual files, we need to convert them to WSG
        # format FRG. For fasta convertTo

        for fasta_fpath, qual_fpath in fasta_qual_files:
            file_info = _scrape_info_from_fname(fasta_fpath)
            library   = file_info['lb']
            platform  = file_info['pl']

            cmd = ['convert-fasta-to-v2.pl', '-noobt', '-l', library, '-s',
                   fasta_fpath, '-q', qual_fpath]

            if platform == '454':
                cmd.append('-454')
            stdout = call(cmd, raise_on_error=True)[0]
            basename  = os.path.splitext(os.path.basename(fasta_fpath))[0]
            frg_fpath = os.path.join(tempdir.name, basename + '.frg')
            frg_fhand = open(frg_fpath, 'w')
            frg_fhand.write(stdout)
            frg_fhand.close()

        # Finally we do a cat of all the frg files ant write it to real output
        # dir
        output_dir = self._create_output_dirs()['assembly_input']
        frg_fhands = []
        for frg in os.listdir(tempdir.name):
            if frg.endswith('.frg'):
                frg_fhands.append(open(os.path.join(tempdir.name, frg)))
        final_frg = os.path.join(output_dir, BACKBONE_BASENAMES['merged_frg'])
        cat(frg_fhands, open(final_frg, 'w'))

        #close all files
        for frg_fhand in frg_fhands:
            frg_fhand.close()
        tempdir.close()

class WSGAssemblyAnalyzer(Analyzer):
    'It assembles the cleaned reads using WSG assembler'

    def run(self):
        '''It runs the analysis.'''
        proj_name = self._get_project_name()
        #we need an assembly dir for this run
        #in this case every assembly is done in a new directory
        #the directory for this assembly
        output_dirs = self._create_output_dirs(timestamped=True)
        assembly_dir = output_dirs['analysis']
        wsg_dir  = os.path.join(assembly_dir, '%s_assembly' % proj_name)

        frgs = self._get_input_fpaths()['frgs']

        stdout = open(os.path.join(assembly_dir, 'stdout.txt'), 'w')
        stderr = open(os.path.join(assembly_dir, 'stderr.txt'), 'w')
        cmd = ['runCA', '-d', wsg_dir, '-p', proj_name]
        cmd.extend(frgs)
        #print ' '.join(cmd)
        raw_input()
        retcode = call(cmd, stdout=stdout, stderr=stderr)[-1]
        if retcode:
            stdout.flush()
            stdout.seek(0)
            raise RuntimeError(stdout.read())
        stdout.close()
        stderr.close()

        # symlinks result to result dir
        fname = '%s.asm' % proj_name
        asm_result_fpath = os.path.join(output_dirs['analysis'], fname)
        os.symlink(asm_result_fpath, os.path.join(output_dirs['result'], fname))

class PrepareMiraAssemblyAnalyzer(Analyzer):
    'It assembles the cleaned reads'

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''

        files_454 = []
        files_sanger_with_qual = []
        files_sanger_without_qual = []
        for fpath in self._get_input_fpaths()['reads']:
            fhand = open(fpath)
            fname = os.path.split(fpath)[-1]
            if 'pl_454' in fname.lower():
                files_454.append(fhand)
            elif 'pl_sanger' in fname.lower():
                format_ = guess_seq_file_format(fhand)
                if format_ == 'fasta':
                    files_sanger_without_qual.append(fhand)
                elif format_ == 'fastq':
                    files_sanger_with_qual.append(fhand)

        #fastq are processed before
        files_sanger = files_sanger_with_qual[:]
        files_sanger.extend(files_sanger_without_qual)

        #all files should be fasta and fasta.qual
        output_dir = self._create_output_dirs()['assembly_input']
        project_name = self._get_project_name()
        for ext, files in (('_in.454', files_454),
                           ('_in.sanger', files_sanger)):
            base_name = os.path.join(output_dir, project_name + ext)
            fasta_fpath = base_name + '.fasta'
            qual_fpath = base_name + '.fasta.qual'
            if os.path.exists(fasta_fpath) or not files:
                continue
            fasta_fhand = open(fasta_fpath, 'w')
            qual_fhand = open(qual_fpath, 'w')
            self._cat_to_fasta(files, fasta_fhand, qual_fhand)
            fasta_fhand.close()
            qual_fhand.close()

        # close all files
        for file_ in files_454 + files_sanger:
            file_.close()

    @staticmethod
    def _files_to_temp_fasta(files):
        'It converts the given files to a temporary fasta and qual'
        fastas, quals = [], []
        for file_ in files:
            #are we dealing with a fastq file (with qual)
            if 'fastq' in os.path.splitext(file_.name)[-1]:
                qual  = NamedTemporaryFile(suffix='.qual')
                fasta = NamedTemporaryFile(suffix='.fasta')
                seqio(in_seq_fhand=file_, out_seq_fhand=fasta,
                        out_qual_fhand=qual, out_format='fasta')
            else:
                #the file is already fasta
                fasta = file_
                qual = None
            fastas.append(fasta)
            quals.append(qual)
        return fastas, quals

    def _cat_to_fasta(self, files, fasta_fhand, qual_fhand):
        'It cats the given files together and it returns a fasta and qual'
        #all files should be in fasta and qual
        fastas, quals = self._files_to_temp_fasta(files)
        cat(fastas, fasta_fhand)
        cat(quals , qual_fhand)

class MiraAssemblyAnalyzer(Analyzer):
    'It assembles the cleaned reads'

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''

        proj_name = self._get_project_name()
        #we need an assembly dir for this run
        #in this case every assembly is done in a new directory
        #the directory for this assembly
        output_dirs = self._create_output_dirs(timestamped=True)
        assembly_dir = output_dirs['analysis']
        mira_dir  = os.path.join(assembly_dir, '%s_assembly' % proj_name)

        os.chdir(assembly_dir)
        #with technologies have we used?
        #each input file should have a link in the assembly dir
        techs = set()
        for fpath in self._get_input_fpaths()['reads']:
            fname = os.path.split(fpath)[-1]
            mira_input_fname = os.path.join(assembly_dir, fname)
            os.symlink(fpath, mira_input_fname)
            if '454' in fname:
                techs.add('454')
            elif 'sanger' in fname:
                techs.add('sanger')

        settings = self._project_settings['Mira']
        #job part of the command

        job = settings['job_options']
        job.extend(techs)
        job = ','.join(job)
        job = '-job=' + job


        cmd = ['mira', '-project=%s' % proj_name, '-fasta', job]
        cmd.append('-OUT:rld=yes')
        for tech in techs:
            tech_str = '%s_settings' % tech
            if tech_str in settings:
                cmd.append(tech_str.upper())
                cmd.extend(settings[tech_str])

        stdout = open(os.path.join(assembly_dir, 'stdout.txt'), 'w')
        stderr = open(os.path.join(assembly_dir, 'stderr.txt'), 'w')
        retcode = call(cmd, stdout=stdout, stderr=stderr)[-1]
        if retcode:
            stdout.flush()
            stdout.seek(0)
            raise RuntimeError(stdout.read())

        stdout.close()
        stderr.close()
        #remove the log directory
        mirachkpt_dir = os.path.join(mira_dir,  "%s_d_chkpt" % proj_name)
        if os.path.exists(mirachkpt_dir):
            shutil.rmtree(mirachkpt_dir)

        #results
        mira_results_dir = os.path.join(mira_dir, '%s_d_results' % proj_name)
        results_dir = output_dirs['result']
        mira_basename = '%s_out' % proj_name
        file_exts = (('.ace', '.ace'),
                     ('.caf', '.caf'),
                     ('.unpadded.fasta', '.fasta'),
                     ('.unpadded.fasta.qual', '.qual'))
        for mira_ext, result_ext in file_exts:
            mira_fpath = os.path.join(mira_results_dir,
                                      mira_basename + mira_ext)
            res_fpath = os.path.join(results_dir,
                                     BACKBONE_BASENAMES['contigs'] + result_ext)
            if os.path.exists(mira_fpath):
                os.symlink(mira_fpath, res_fpath)


        mira_info_dir = os.path.join(mira_dir, '%s_d_info' % proj_name)
        if os.path.exists(mira_info_dir):
            os.symlink(mira_info_dir, os.path.join(results_dir,
                                                  BACKBONE_DIRECTORIES['info']))

class LastAnalysisAnalyzer(Analyzer):
    'It chooses the latest assembly as the result'

    def run(self):
        '''It runs the analysis. It checks if the analysis'''
        return self._select_last_analysis()

    def _select_last_analysis(self):
        '''It select the last analysis from a directory with timestamped results'''
        dirs = self._get_timestamped_output_dirs()
        latest_dir = dirs[-1]
        output_dir = self._get_output_dirs()['result']
        if os.path.exists(output_dir):
            os.remove(output_dir)
        os.symlink(latest_dir, output_dir)

class SetAssemblyAsReferenceAnalyzer(Analyzer):
    'It sets the reference assembly as mapping reference'
    def run(self):
        '''It runs the analysis.'''
        contigs_fpath = self._get_input_fpaths()['contigs']
        contigs_ext = os.path.splitext(contigs_fpath)[-1]
        reference_dir = self._create_output_dirs()['result']
        reference_fpath = os.path.join(reference_dir,
                          BACKBONE_BASENAMES['mapping_reference'] + contigs_ext)
        os.symlink(contigs_fpath, reference_fpath)

def _get_basename(fpath):
    'It returns the base name without path and extension'
    return os.path.splitext(os.path.basename(fpath))[0]

class MappingAnalyzer(Analyzer):
    'It performs the mapping of the sequences to the reference'
    def run(self):
        '''It runs the analysis.'''
        settings = self._project_settings['Mappers']
        inputs = self._get_input_fpaths()
        reads_fpaths = inputs['reads']
        reference_fpath = inputs['reference']
        output_dir = self._create_output_dirs(timestamped=True)['result']

        for read_fpath in reads_fpaths:
            read_info = _scrape_info_from_fname(read_fpath)
            platform = read_info['pl']
            #which mapper are we using for this platform
            mapper = settings['mapper_for_%s' % platform]
            out_bam = os.path.join(output_dir,
                                   _get_basename(read_fpath) + '.bam')
            mapping_parameters = {}
            if platform in ('454', 'sanger'):
                mapping_parameters['reads_length'] = 'long'
            else:
                mapping_parameters['reads_length'] = 'short'
            map_reads(mapper,
                      reads_fpath=read_fpath,
                      reference_fpath=reference_fpath,
                      out_bam_fpath=out_bam,
                      parameters = mapping_parameters)

class MergeBamAnalyzer(Analyzer):
    'It performs the merge of various bams into only one'
    def run(self):
        '''It runs the analysis.'''
        settings = self._project_settings
        project_path = settings['General_settings']['project_path']
        os.chdir(project_path)
        inputs          = self._get_input_fpaths()
        bam_fpaths      = inputs['bams']
        reference_fpath = inputs['reference']

        output_dir      = self._create_output_dirs()['result']
        merged_bam_fpath = os.path.join(output_dir,
                                       BACKBONE_BASENAMES['merged_bam'])
        # First we need to create the sam with added tags and headers
        temp_dir = NamedTemporaryDir()
        for bam_fpath in bam_fpaths:
            bam_basename = os.path.splitext(os.path.basename(bam_fpath))[0]
            temp_sam     =  NamedTemporaryFile(prefix='%s.' % bam_basename,
                                               suffix='.sam')
            sam_fpath    = os.path.join(temp_dir.name, bam_basename + '.sam')
            bam2sam(bam_fpath, temp_sam.name)
            sam_fhand = open(sam_fpath, 'w')
            add_header_and_tags_to_sam(temp_sam, sam_fhand)
            sam_fhand.close()

        # Once the headers are ready we are going to merge
        sams = []
        for file_ in os.listdir(temp_dir.name):
            if file_.endswith('.sam'):
                sams.append(open(os.path.join(temp_dir.name, file_)))

        temp_sam = NamedTemporaryFile(suffix='.sam')
        reference_fhand = open(reference_fpath)
        merge_sam(sams, temp_sam, reference_fhand)
        reference_fhand.close()

        # close files
        for sam in sams:
            sam.close()
        # Convert sam into a bam,(Temporary)
        temp_bam = NamedTemporaryFile(suffix='.bam')
        sam2bam(temp_sam.name, temp_bam.name)

        # finally we need to order the bam
        sort_bam_sam(temp_bam.name, merged_bam_fpath)
        temp_bam.close()
        temp_sam.close()

def _get_basename(fpath):
    'It returns the basename for the given path'
    return os.path.splitext(os.path.basename(fpath))[0]

def _get_seq_or_repr_fpath(seqs_fpaths, repr_fpaths):
    'It returns for every file the repr or the seq file'
    repr_fpaths = dict([(_get_basename(fpath), fpath) for fpath in repr_fpaths])
    new_seq_fpaths = []
    for fpath in seqs_fpaths:
        basename = _get_basename(fpath)
        if basename in repr_fpaths:
            new_seq_fpaths.append(repr_fpaths[basename])
        else:
            new_seq_fpaths.append(fpath)
    return new_seq_fpaths

class SnvCallerAnalyzer(Analyzer):
    'It performs the calling of the snvs in a bam file'

    def run(self):
        'It runs the analysis.'
        output_dir   = self._create_output_dirs()['result']
        inputs       = self._get_input_fpaths()
        repr_fpaths  = inputs['repr']
        seqs_fpaths  = inputs['input']
        seqs_fpaths  = _get_seq_or_repr_fpath(seqs_fpaths, repr_fpaths)
        merged_bam   = inputs['merged_bam']

        pipeline = 'snv_bam_annotator'
        bam_fhand = open(merged_bam)
        configuration = {'snv_bam_annotator': {'bam_fhand':bam_fhand}}
        settings = self._project_settings
        if 'Snvs' in settings and 'min_quality' in settings['Snvs']:
            min_quality = settings['Snvs']['min_quality']
            configuration['snv_bam_annotator']['min_quality'] = int(min_quality)

        for seq_fpath in seqs_fpaths:
            temp_repr = NamedTemporaryFile(suffix='.repr', mode='a',
                                           delete=False)
            io_fhands = {'in_seq': open(seq_fpath),
                         'outputs':{'repr':temp_repr}}
            seq_pipeline_runner(pipeline, configuration=configuration,
                                io_fhands=io_fhands)
            temp_repr.close()
            repr_fpath = os.path.join(output_dir,
                                      _get_basename(seq_fpath) + '.repr')
            shutil.move(temp_repr.name, repr_fpath)

class WriteAnnotationAnalyzer(Analyzer):
    'It writes all the output annoation files'
    def run(self):
        'It runs the analysis.'
        output_dir   = self._create_output_dirs()['result']
        inputs       = self._get_input_fpaths()
        repr_fpaths  = inputs['repr']

        output_files = ['vcf']
        for seq_fpath in repr_fpaths:
            outputs = {}
            for output_kind in output_files:
                output_fpath = os.path.join(output_dir,
                                   _get_basename(seq_fpath) + '.' + output_kind)
                if os.path.exists(output_fpath):
                    os.remove(output_fpath)
                output_fhand = open(output_fpath, 'a')
                outputs[output_kind] = output_fhand
            io_fhands = {'in_seq': open(seq_fpath),
                         'outputs': outputs}
            seq_pipeline_runner(pipeline=None,
                                configuration=None,
                                io_fhands=io_fhands)

ANALYSIS_DEFINITIONS = {
    'clean_reads':
        {'inputs':{
            'reads':
                {'directory': 'original_reads',
                 'file_kinds': 'sequence_files'}
            },
         'outputs':{'reads':{'directory': 'cleaned_reads'}},
         'analyzer': CleanReadsAnalyzer,
        },
    'prepare_wsg_assembly':
        {'inputs':{
            'reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'}
            },
         'outputs':{'assembly_input':{'directory': 'assembly_input'}},
         'analyzer': PrepareWSGAssemblyAnalyzer,
        },
    'wsg_assembly':
        {'inputs':{
            'frgs':
                {'directory': 'assembly_input',
                 'file_kinds': 'frg_file'}
            },
         'outputs':{
                    'analysis':  {'directory': 'assemblies'},
                    'result':  {'directory': 'assembly_result'},
                    },
         'analyzer': WSGAssemblyAnalyzer,
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
    'prepare_mira_assembly':
        {'inputs':{
            'reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'}
            },
         'outputs':{'assembly_input':{'directory': 'assembly_input'}},
         'analyzer': PrepareMiraAssemblyAnalyzer,
        },
    'mira_assembly':
        {'inputs':{
            'reads':
                {'directory': 'assembly_input',
                 'file_kinds': 'sequence_or_qual_files'}
            },
         'outputs':{
                    'analysis':  {'directory': 'assemblies'},
                    'result':  {'directory': 'assembly_result'},
                    },
         'analyzer': MiraAssemblyAnalyzer,
        },
    'select_last_assembly':
        {'inputs':{'analyses_dir':{'directory': 'assemblies'}},
         'outputs':{'result':{'directory': 'assembly_result',
                              'create': False}},
         'analyzer': LastAnalysisAnalyzer,
        },
    'set_assembly_as_reference':
        {'inputs':{
                   'contigs':
                            {'directory': 'assembly_result',
                             'file': 'contigs'},
                   },
         'outputs':{'result':{'directory': 'mapping_reference'}},
         'analyzer': SetAssemblyAsReferenceAnalyzer,
        },
    'mapping':
        {'inputs':{
            'reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'},
            'reference':
                {'directory': 'mapping_reference',
                'file': 'mapping_reference'},

            },
         'outputs':{'result':{'directory': 'mappings_by_readgroup'}},
         'analyzer': MappingAnalyzer,
        },
    'select_last_mapping':
        {'inputs':{'analyses_dir':{'directory': 'mappings'}},
         'outputs':{'result':{'directory': 'mapping_result',
                              'create':False}},
         'analyzer': LastAnalysisAnalyzer,
        },
    'merge_bam':
        {'inputs':{
            'bams':
                {'directory': 'mappings_by_readgroup',
                 'file_kinds': 'bam'},
            'reference':
                {'directory': 'mapping_reference',
                'file': 'mapping_reference'},
            },
         'outputs':{'result':{'directory': 'mapping_result'}},
         'analyzer': MergeBamAnalyzer,
        },
    'call_snv':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'merged_bam':
                {'directory':'mapping_result',
                 'file': 'merged_bam'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': SnvCallerAnalyzer},
    'write_annotation':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_result'}},
         'analyzer': WriteAnnotationAnalyzer},
}

BACKBONE_DIRECTORIES = {
    'config_file': 'backbone.conf',
    'original_reads': 'reads/original',
    'cleaned_reads': 'reads/cleaned',
    'assembly_input': 'assembly/input',
    'assemblies': ('assembly', ''),
    'assembly_result': ('assembly', 'result'),
    'mappings': ('mapping', ''),
    'mapping_result': ('mapping', 'result'),
    'mapping_reference': 'mapping/reference',
    'mappings_by_readgroup': ('mapping', 'result/by_readgroup'),
    'pileups':'mapping/result/pileups',
    'snvs':'annotations/snvs',
    'info':'info',
    'original_reads_stats': 'reads/original/stats',
    'cleaned_reads_stats': 'reads/cleaned/stats',
    'annotation_repr':'annotations/repr',
    'annotation_input':'annotations/input',
    'annotation_result':'annotations/result',
                       }
BACKBONE_BASENAMES = {
    'contigs':'contigs',
    'mapping_reference':'reference',
    'merged_bam':'merged.bam',
    'snv_result':'all.snvs',
    'merged_frg':'all_seq.frg',
}

def do_analysis(kind, project_settings=None, analysis_config=None):
    'It does one of the predefined analyses'
    if project_settings is None:
        project_settings = os.path.join(os.getcwd(),
                                        BACKBONE_DIRECTORIES['config_file'])
        if not os.path.exists(project_settings):
            raise ValueError('Settings path not given and not found')

    if not analysis_config:
        analysis_config = {}

    settings = ConfigObj(project_settings)
    analysis_def = ANALYSIS_DEFINITIONS[kind]

    analyzer_klass = analysis_def['analyzer']
    analyzer = analyzer_klass(project_settings=settings,
                        analysis_definition=analysis_def)

    analyzer.run()
