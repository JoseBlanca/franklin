'''
This module is part of ngs_backbone. This is the general analysis repository.
They do not perform any analysis but proveds of genearl methosd for other
analysis.

Created on 01/03/2010

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

import os, time, tempfile, logging
import franklin
from franklin.seq.readers import guess_seq_file_format
from franklin.backbone.specifications import (BACKBONE_DIRECTORIES,
                                              BACKBONE_BASENAMES)
from franklin.utils.misc_utils import VersionedPath

def scrape_info_from_fname(path):
    'It guess pipeline taking into account the platform and the file format'
    file_info = {}
    fhand = open(path.last_version)
    file_info['format'] = guess_seq_file_format(fhand)
    fhand.close()
    fname = path.basename
    for item in fname.split('.'):
        key, value = item.split('_', 1)
        file_info[key] = value
    return file_info

def _is_sequence_file(path):
    'It returns true if the function is a sequence'
    return path.extension in ('fasta', 'fastq', 'sfastq', 'repr')

def _is_sequence_or_qual_file(path):
    'It returns true if the function is a sequence or quality file'
    return path.extension in ('fasta', 'fastq', 'sfastq', 'repr', 'qual')

def _select_fname(kind, paths):
    'It returns the path that correponds to the given file kind'
    if not paths:
        raise RuntimeError('No paths given, so no selection possible')

    conf_fname = BACKBONE_BASENAMES[kind]
    conf_basename, conf_ext = os.path.splitext(conf_fname)
    conf_ext = conf_ext[1:]
    if kind in ('contigs', 'mapping_reference'):
        paths = filter(_is_sequence_file, paths)
        paths = filter(lambda x: conf_basename == x.basename, paths)
    elif kind in ('merged_bam'):
        paths = filter(lambda x: conf_basename == x.basename, paths)
        paths = filter(lambda x: conf_ext == x.extension, paths)

    if not paths:
        raise RuntimeError('No file found for %s' % kind)
    assert len(paths) == 1
    return paths[0]

class Analyzer(object):
    '''This class performs an analysis.

    This is a prototype that should be inherited for each analysis. This class
    looks for inputs files and runs the analysis.
    It checks if the the analysis is already done, looking to the output
    directory
    '''
    def __init__(self, project_settings, analysis_definition, silent=False):
        'The init'
        self._project_settings = project_settings
        self._analysis_def = analysis_definition
        self._timestamped_dir = None
        self._old_tmpdir = tempfile.gettempdir()
        self._setup_tempdir()
        self._silent = silent

        if 'threads' in self._project_settings['General_settings']:
            self.threads = self._project_settings['General_settings']['threads']
        else:
            self.threads = False

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
            timed_dir = os.path.abspath(timed_dir)
            dirs[kind] = timed_dir
        return dirs

    def _get_input_fpaths(self):
        'It returns a dict with the inputs files'
        inputs_def = self._analysis_def['inputs']
        input_files = {}
        input_dirs = self._get_input_dirs()
        for input_kind, input_def in inputs_def.items():
            backbone_dir = input_dirs[input_kind]
            paths = VersionedPath(backbone_dir).list_paths_versioned()

            if 'file_kinds' in input_def:
                if input_def['file_kinds'] == 'sequence_files':
                    paths = filter(_is_sequence_file, paths)
                elif input_def['file_kinds'] == 'sequence_or_qual_files':
                    paths = filter(_is_sequence_or_qual_file, paths)
                elif input_def['file_kinds'] == 'bam':
                    paths = filter(lambda x: x.extension == 'bam', paths)
                elif input_def['file_kinds'] == 'frg_file':
                    paths = filter(lambda x: x.extension == 'frg', paths)
            if 'file' in input_def:
                paths = _select_fname(input_def['file'], paths)
            input_files[input_kind] = paths
        return input_files

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

    _analysis_messages = {
        'AnnotateIntronsAnalyzer': 'Annotating introns',
        'AnnotateDescriptionAnalyzer': 'Annotating descriptions',
        'AnnotateGoAnalyze': 'Annotating GO terms',
        'AnnotateMicrosatelliteAnalyzer': 'Annotating Microsatellites',
        'AnnotateOrfAnalyzer': 'Annotating ORFs',
        'AnnotateOrthologsAnalyzer': 'Annotating Orthologs',
                          }

    def _log(self, messages):
        'It logs the analysis'
        if not self._silent:
            #create logger
            logger = logging.getLogger("franklin")
            #which analysis is running?
            class_name = self.__class__.__name__
            class_name = class_name.split('.')[-1]
            if class_name in self._analysis_messages:
                analysis_message = self._analysis_messages[class_name]
            else:
                analysis_message = class_name
            if 'analysis_started' in messages:
                logger.info(analysis_message)
                logger.info('backbone version: %s' % str(franklin.__version__))
                logger.info('Analysis started')
            elif 'analysis_finished' in messages:
                logger.info('Analysis finished')

class LastAnalysisAnalyzer(Analyzer):
    'It chooses the latest assembly as the result'

    def run(self):
        '''It runs the analysis. It checks if the analysis'''
        return self._select_last_analysis()

    def _select_last_analysis(self):
        '''It select the last analysis from a directory with timestamped results'''
        self._log({'analysis_started':True})
        dirs = self._get_timestamped_output_dirs()
        latest_dir = dirs[-1]
        output_dir = self._get_output_dirs()['result']
        if os.path.exists(output_dir):
            os.remove(output_dir)
        os.symlink(latest_dir, output_dir)
