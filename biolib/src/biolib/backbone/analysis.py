'''
Created on 29/01/2010

@author: jose
'''

import os, tempfile, time
from biolib.utils.seqio_utils import guess_seq_file_format, seqio, cat
from biolib.pipelines.pipelines import seq_pipeline_runner
from biolib.utils.cmd_utils import call
from configobj import ConfigObj

def _is_sequence_file(fpath):
    'It returns true if the function is a sequence'
    ext = os.path.splitext(fpath)[-1].strip('.')
    if ext in ('fasta', 'fastq', 'sfastq'):
        return True
    else:
        return False


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
        self._create_output_dir()

    def _get_input_fpaths(self):
        'It returns a dict with the inputs files'
        inputs_def = self._analysis_def['inputs']
        input_files = {}
        for input_kind, input_def in inputs_def.items():
            project_dir = self._project_settings['General_settings']['project_path']
            backbone_dir = os.path.join(project_dir,
                                   BACKBONE_DIRECTORIES[input_def['directory']])
            fpaths = os.listdir(backbone_dir)

            #we want full paths
            fpaths = [os.path.join(backbone_dir, path) for path in fpaths]

            if input_def['file_kinds'] == 'sequence_files':
                fpaths = filter(_is_sequence_file, fpaths)
            input_files[input_kind] = fpaths
        return input_files

    def get_output_fpaths(self):
        'It returns the output files for the analysis'
        raise NotImplementedError

    def run(self):
        'It runs the analysis'
        raise NotImplementedError

    def _get_output_dir(self):
        'It returns the path for the output dir'
        dir_kind = self._analysis_def['output']['directory']
        project_dir = self._project_settings['General_settings']['project_path']
        output_dir = os.path.join(project_dir,
                                  BACKBONE_DIRECTORIES[dir_kind])
        return output_dir

    def _create_output_dir(self):
        'It creates the output directory if it does not exists'
        output_dir = self._get_output_dir()
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

class CleanReadsAnalyzer(Analyzer):
    'It does a cleaning reads analysis'
    def _get_output_fpath(self, input_fpath):
        'It returns the output files for the analysis'
        ouput_dir_kind = self._analysis_def['output']['directory']
        project_path = self._project_settings['General_settings']['project_path']
        output_dir = os.path.join(project_path,
                                  BACKBONE_DIRECTORIES[ouput_dir_kind])
        fname = os.path.split(input_fpath)[-1]
        return  os.path.join(output_dir, fname)

    def _guess_cleaning_pipepile(self, file_info):
        'It returns the pipeline suited to clean the given file'
        if file_info['format'] == 'fasta':
            return 'sanger_without_qual'
        elif file_info['pt'] == 'illumina':
            return 'solexa'
        elif file_info['pt'] == '454':
            return 'sanger_with_qual'
        elif file_info['pt'] == 'sanger':
            return 'sanger_with_qual'
        else:
            raise ValueError('Unable to guess the cleaning pipeline: %s' %
                             str(file_info))

    def create_cleaning_configuration(self, platform):
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

        #TODO lucy settings
        return configuration

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        inputs   = self._get_input_fpaths()['reads']

        for input_fpath in inputs:
            output_fpath = self._get_output_fpath(input_fpath)
            if os.path.exists(output_fpath):
                continue
            file_info = _scrape_info_from_fname(input_fpath)
            pipeline = self._guess_cleaning_pipepile(file_info)
            iofhands = {'in_seq':open(input_fpath),
                        'out_seq': open(output_fpath, 'w')}
            configuration = self.create_cleaning_configuration(
                                                       platform=file_info['pt'])
            seq_pipeline_runner(pipeline, configuration, iofhands,
                                file_info['format'])
        return

def _scrape_info_from_fname(fpath):
    'It guess pipeline taking into account the platform and the file format'
    file_info = {}
    file_info['format'] = guess_seq_file_format(open(fpath, 'r'))

    fname = os.path.splitext(os.path.split(fpath)[-1])[0]
    for item in fname.split('.'):
        key, value     = item.split('_', 1)
        file_info[key] = value
    return file_info

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
            if 'pt_454' in fname.lower():
                files_454.append(fhand)
            elif 'pt_sanger' in fname.lower():
                format_ = guess_seq_file_format(fhand)
                if format_ == 'fasta':
                    files_sanger_without_qual.append(fhand)
                elif format_ == 'fastq':
                    files_sanger_with_qual.append(fhand)
        #fastq are processed before
        files_sanger = files_sanger_with_qual[:]
        files_sanger.extend(files_sanger_without_qual)

        #all files should be fasta and fasta.qual
        for ext, files in (('_in.454', files_454), ('_in.sanger', files_sanger)):
            base_name = os.path.join(self._get_output_dir(),
                                   self._project_settings['General_settings']['project_name'] + ext)
            fasta_fpath = base_name + '.fasta'
            qual_fpath = base_name + '.fasta.qual'
            if os.path.exists(fasta_fpath) or not files:
                continue
            fasta_fhand = open(fasta_fpath, 'w')
            qual_fhand = open(qual_fpath, 'w')
            self._cat_to_fasta(files, fasta_fhand, qual_fhand)

    @staticmethod
    def _files_to_temp_fasta(files):
        'It converts the given files to a temporary fasta and qual'
        fastas, quals = [], []
        for file_ in files:
            #are we dealing with a fastq file (with qual)
            if 'fastq' in os.path.splitext(file_.name)[-1]:
                qual = tempfile.NamedTemporaryFile(suffix='.qual')
                fasta = tempfile.NamedTemporaryFile(suffix='.fasta')
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

    def _create_assembly_dir(self, result_dir):
        'It returns a directory for this analysis'
        assembly_dir = os.path.split(result_dir.rstrip('/'))[0]
        assembly_dir = os.path.join(assembly_dir,
                                    time.strftime("%Y%m%d_%H%M", time.gmtime()))
        os.mkdir(assembly_dir)
        return assembly_dir

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        output_dir = self._get_output_dir()

        #we need an assembly dir for this run
        #in this case every assembly is done in a new directory
        assembly_dir = self._create_assembly_dir(output_dir)

        #the directory for this assembly

        os.chdir(assembly_dir)
        #with technologies have we used?
        #each input file should have a link in the assembly dir
        techs = set()
        for fpath in self._get_input_fpaths()['reads']:
            fname = os.path.split(fpath)[-1]
            mira_input_fname = os.path.join(assembly_dir, fname)
            os.symlink(fname, mira_input_fname)
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
        project_name = self._project_settings['General_settings']['project_name']

        cmd = ['mira', '-project=%s' % project_name, '-fasta', job]

        for tech in techs:
            tech_str = '%s_settings' % tech
            if tech_str in settings:
                cmd.append(tech_str.upper())
                cmd.extend(settings[tech_str])

        stdout = open(os.path.join(assembly_dir, 'stdout.txt'), 'w')
        stderr = open(os.path.join(assembly_dir, 'stderr.txt'), 'w')
        retcode = call(cmd, stdout=stdout, stderr=stderr)[-1]

        #remove the log directory

        #link the results to the result directory
        return retcode



ANALYSIS_DEFINITIONS = {
    'clean_reads':
        {'inputs':{
            'reads':
                {'directory': 'original_reads',
                 'file_kinds': 'sequence_files'}
            },
         'output':{'directory': 'cleaned_reads'},
         'analyzer': CleanReadsAnalyzer,
        },
    'prepare_mira_assembly':
        {'inputs':{
            'reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'}
            },
         'output':{'directory': 'assembly_input'},
         'analyzer': PrepareMiraAssemblyAnalyzer,
        },
    'mira_assembly':
        {'inputs':{
            'reads':
                {'directory': 'assembly_input',
                 'file_kinds': 'sequence_files'}
            },
         'output':{'directory': 'assembly_result'},
         'analyzer': MiraAssemblyAnalyzer,
        },
}

BACKBONE_DIRECTORIES = {
    'config_file': 'backbone.conf',
    'original_reads': 'reads/original',
    'cleaned_reads': 'reads/cleaned',
    'assembly_input': 'assembly/input',
    'assembly_result': 'assembly/result',
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
