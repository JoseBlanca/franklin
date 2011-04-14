'''
Command utilities for franklin

This module provides utilities to run external commands into franklin
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

import subprocess, signal, tempfile, os, itertools
import StringIO, logging, copy, shutil, platform

from franklin.seq.writers import temp_fasta_file, temp_qual_file
from franklin.utils.misc_utils import (NamedTemporaryDir, DisposableFile,
                                       get_franklin_ext_dir)

def _locate_file(fpath):
    cmd = ['locate', fpath]
    stdout, stderr, retcode = call(cmd, raise_on_error=False, add_ext_dir=False)
    if retcode == 1:
        raise RuntimeError('File not found: %s' % fpath)
    elif not retcode:
        pass
    else:
        msg = 'Locate produced and error when looking for: %s' % fpath
        raise RuntimeError(msg)
    found_path = None
    for line in stdout.splitlines():
        if fpath in line:
            found_path = line.strip()
            break
    return found_path

def _guess_java_install_dir(jar_fpath):
    'It returns the dir path using locate on a jar file'
    java_dir_path = _locate_file(jar_fpath)
    if not java_dir_path:
        msg = '%s was not found in your system and it is required' % jar_fpath
        raise RuntimeError(msg)
    java_dir_path = java_dir_path.replace(jar_fpath, '')
    return java_dir_path

def guess_jar_dir(jar_name, java_conf=None):
    'It returns the jar_name path using locate'
    conf_variables = {'SortSam.jar': 'picard_path',
                      'GenomeAnalysisTK.jar': 'gatk_path',
                      'blast2go.jar': 'blast2go_path'}
    java_dir_names = {'SortSam.jar': 'picard',
                      'GenomeAnalysisTK.jar': 'gatk',
                      'blast2go.jar': 'blast2go'}
    jar_path = None
    if jar_name in conf_variables:
        conf_var = conf_variables[jar_name]
        if java_conf and conf_var in java_conf and java_conf[conf_var]:
            jar_path = java_conf[conf_var]
    if not jar_path:
        if jar_name == 'blast2go.jar':
            jar_path = _guess_java_install_dir(jar_name)
        else:
            franklin_path = get_franklin_ext_dir()
            jar_path = os.path.join(franklin_path, 'java',
                                    java_dir_names[jar_name])

    return jar_path


# Runner definitions, Define here the parameters of the programs you want to
# use with this class
STDOUT = 'stdout'
ARGUMENT = 'argument'
STDIN = 'stdin'

BLASTPLUS_DEF = {'binary':'',
                 'parameters': {'database' :{'option': '-db'},
                   'expect'   :   {'default': 0.0001,'option': '-evalue'},
                   'nhitsv'   :   {'default': 20, 'option':'-num_descriptions'},
                   'nhitsb'   :   {'default': 20, 'option':'-num_alignments'},
                   'alig_format': {'default':5, 'option':'-outfmt'},
                   'gapextend':   {'option': '-gapextend'},
                   'gapopen':     {'option': '-gapopen'},
                   'task':        {'option': '-task'},
                   'subject':     {'option': '-subject'},
                   'no_greedy':   {'option': '-no_greedy'}
                            },
                 'output':{'blast+':{'option':STDOUT}},
            'input':{'sequence':{'option':'-query', 'files_format':['fasta']}},
                 'ignore_stderrs': ['Karlin-Altschul']}

BLASTN_DEF = copy.deepcopy(BLASTPLUS_DEF)
BLASTN_DEF['binary'] = 'blastn'
BLASTN_DEF['output'] = {'blastn':{'option':STDOUT}}
BLASTN_DEF['parameters']['penalty'] = {'option':'-penalty'}
BLASTN_DEF['parameters']['dust'] = {'option':'-dust'}
BLASTN_DEF['parameters']['reward'] = {'option':'-reward'}
BLASTN_DEF['parameters']['penalty'] = {'option':'-penalty'}

BLASTP_DEF = copy.deepcopy(BLASTPLUS_DEF)
BLASTP_DEF['binary'] = 'blastp'
BLASTP_DEF['output'] = {'blastp':{'option':STDOUT}}

BLASTX_DEF = copy.deepcopy(BLASTPLUS_DEF)
BLASTX_DEF['binary'] = 'blastx'
BLASTX_DEF['output'] = {'blastx':{'option':STDOUT}}

TBLASTN_DEF = copy.deepcopy(BLASTPLUS_DEF)
TBLASTN_DEF['binary'] = 'tblastn'
TBLASTN_DEF['output'] = {'tblastn':{'option':STDOUT}}

TBLASTX_DEF = copy.deepcopy(BLASTPLUS_DEF)
TBLASTX_DEF['binary'] = 'tblastx'
TBLASTX_DEF['output'] = {'tblastx':{'option':STDOUT}}

#TODO megablast?


RUNNER_DEFINITIONS = {
    'blastn':BLASTN_DEF,
    'blastp':BLASTP_DEF,
    'blastx':BLASTX_DEF,
    'tblastn':TBLASTN_DEF,
    'tblastx':TBLASTX_DEF,
    'water': {'binary':'water',
            'parameters': {'subject' :{'required':True, 'option': '-bsequence'},
                           'gapopen'   :{'default': 10.0,'option': '-gapopen'},
                           'gapextend' :{'default': 0.5, 'option':'-gapextend'},
                    'outformat'   :{'default': 'markx10', 'option':'-aformat3'},
                    },
            'output':{'water':{'option':'-outfile'}},
            'input':{'sequence':{'option':'-asequence', 'files_format':['fasta']}}
              },
    'mdust':{'binary':'mdust',
             'parameters':{'mask_letter':{'default':'L', 'option' : '-m'},
                           'cut_off'    :{'default':'25', 'option':'-v' },
                           'show_masked_segments':{'default':None,
                                                  'option':'-c'}},
             'output':{'sequence':{'option':STDOUT}},
             'input':{'sequence':{'option':ARGUMENT, 'arg_before_params':True,
                                  'files_format':['fasta']}}
            },
    'trimpoly':{'binary':'trimpoly',
                'parameters':{'min_score':{'option':'-s'},
                              'end':{'option':'-e'},
                              'incremental_dist':{'option':'-l'},
                              'fixed_dist':{'option':'-L'},
                              'only_n_trim':{'option':'-N'},
                              'ntrim_above_percent':{'option':'-n'}
                              },
             'output':{'sequence':{'option':STDOUT}},
             'input':{'sequence': {'option':STDIN, 'files_format':['fasta']}}
                },
    'lucy':{'binary':'lucy',
            'parameters':{
                      'cdna'   :{'option':'-c', 'default':None},
                      'keep'   :{'option':'-k', 'default':None},
                      'bracket':{'option':'-b', 'default':[10, 0.02]},
                      'window' :{'option':'-w', 'default':[50, 0.08, 10, 0.3]},
                      'error'  :{'option':'-e', 'default':[0.015, 0.015]},
                      'vector' :{'option':'-vector'}
                      },
            'input':{'sequence':{'option': ARGUMENT,
                                 'arg_before_params':True,
                                 'files_format':['fasta', 'qual']}},
            'output':{'sequence':{'option': '-output',
                                  'files_format':['fasta', 'qual']}}
            },
    'sputnik':{'binary':'sputnik',
               'parameters':{
                             'max_unit_length':{'option':'-u', 'default':4},
                             'min_unit_length':{'option':'-v', 'default':2},
                             'min_length_of_ssr':{'option':'-L', 'default':20}
                             },
               'input':{'sequence':{'option':ARGUMENT,
                                    'arg_before_params':False,
                                    'files_format':['fasta']}
                       },
               'output':{'sputnik':{'option':STDOUT}}
               },
    'estscan':{'binary':'estscan',
               'parameters':{ 'matrix':{'option':'-M'}},
               'input':{'sequence':{'option':ARGUMENT,
                                    'arg_before_params':False,
                                    'files_format':['fasta']}},
               'output':{'dna':{'option':STDOUT, 'files_format':['fasta']},
                         'protein':{'option': '-t', 'files_format':['fasta']}}
               },
    }

def _process_parameters(parameters, parameters_def):
    '''Given the parameters definition and some parameters it process the params
    It returns the parameters need by programa'''
    #we process all the parameters
    for param, definition in parameters_def.iteritems():
        #the requiered parameters
        if 'required' in definition and param not in parameters:
            msg = 'parameter ' + param + 'should be given for the cmd'
            raise ValueError(msg)
        #the default parameters
        if 'default' in definition and param not in parameters:
            parameters[param] = definition['default']

    #create the bin for the cmd
    bin_ = []
    if 'bin' in parameters:
        bin_ = parameters['bin']

    for param, value in parameters.items():
        if param == 'bin':
            continue
        param_opt = parameters_def[param]['option']
        bin_.append(param_opt)
        # Values can be a list of parameters
        if isinstance(value, list) or isinstance(value, tuple):
            bin_.extend([ _param_to_str(value_) for value_ in value])
        else:
            if value is not None:
                bin_.append(_param_to_str(value))
    return bin_

def _param_to_str(param):
    'given a parameter It returns an str, that can be used by an CLI program'
    try:
        # If it is a file...
        param = param.name
    except AttributeError:
        pass
    return str(param)

def _prepare_input_files(inputs, seqs):
    'It prepares inputs taking into account the format'
    for key, value in inputs.items():
        files_format = value['files_format']
        inputs[key]['fhands'] = []
        inputs[key]['fpaths'] = []
        seqs, seqs_qual = itertools.tee(seqs, 2)
        for file_format in files_format:
            if file_format == 'fasta':
                fhand = temp_fasta_file(seqs=seqs)
            elif file_format == 'qual':
                fhand = temp_qual_file(seqs=seqs_qual)
            inputs[key]['fhands'].append(fhand)
            inputs[key]['fpaths'].append(fhand.name)

def _get_mktemp_fpaths(num_fpaths):
    'It returns the name of some temp file'
    fpaths = []
    for index in range(num_fpaths):
        fhand, fpath = tempfile.mkstemp()
        # We use here mkstemp because we only need the path, but in order to
        # not have too much open files we have to manualy close the fhand
        os.close(fhand)
        fpaths.append(fpath)
    for fpath in fpaths:
        os.remove(fpath)
    return fpaths

def _prepare_output_files(outputs):
    'It prepares inputs taking into account the format'
    for key, value in outputs.items():
        if 'files_format' in value:
            fpaths = _get_mktemp_fpaths(len(value['files_format']))
        else:
            fpaths = _get_mktemp_fpaths(1)
        outputs[key]['fpaths'] = fpaths

def _build_cmd(cmd_params, runner_def):
    'It bulds the cmd line using  the command definitions'
    inputs = runner_def['input']
    outputs = runner_def['output']
    stdin = None
    bin = runner_def['binary']
    cmd_args_begin = []
    cmd_args_end = []

    for parameters in (inputs, outputs):
        for parameter in parameters.values():
            fpaths = parameter['fpaths']
            if parameter['option'] == STDIN:
                parameter['fhands'][0].seek(0)
                stdin = parameter['fhands'][0].read()
            elif parameter['option'] == STDOUT:
                pass
            elif (parameter['option'] == ARGUMENT and
                  parameter['arg_before_params']):
                cmd_args_begin.extend(fpaths)
            elif (parameter['option'] == ARGUMENT and
                  not parameter['arg_before_params']):
                cmd_args_end.extend(fpaths)
            else:
                cmd_params.append(parameter['option'])
                #know we need to append the output_files
                cmd_params.extend(fpaths)

    cmd = [bin]
    cmd.extend(cmd_args_begin)
    cmd.extend(cmd_params)
    cmd.extend(cmd_args_end)
    return cmd, stdin

def create_runner(tool, parameters=None, environment=None):
    ''''It creates a runner class.

    The runner will be able to run a binary program for different sequences.
    kind is the type of runner (blast, seqclean, etc)
    if multiseq is True the runner will expect list or iterator of sequences.
    '''
    # process parameters to build the cmd
    if parameters is None:
        parameters = {}
    if environment is None:
        environment = {}

    # You can pass the runner definition to create_runner directly in the tool
    # argument
    if isinstance(tool, str):
        runner_def = RUNNER_DEFINITIONS[tool]
    else:
        runner_def = tool
        tool       =  runner_def['binary']
    general_cmd_param = _process_parameters(parameters, runner_def['parameters'])

    def run_cmd_for_sequence(sequence):
        'It returns a result for the given sequence or sequences'
        #parameters should be in the scope because some tempfile could be in
        #there. In some pythons this has been a problem.
        runner_data = copy.deepcopy(runner_def)
        cmd_param = copy.deepcopy(general_cmd_param)
        assert type(parameters)

        #is this a sequence or a generator with seqs
        methods = dir(sequence)
        if 'annotations' in methods or 'lower' in methods:
            sequences = (sequence,)
        else:
            sequences = sequence

        _prepare_input_files(runner_data['input'], sequences)
        _prepare_output_files(runner_data['output'])
        cmd, stdin = _build_cmd(cmd_param, runner_data)

        #print ' '.join(cmd)
        #raw_input()
        try:
            stdout, stderr, retcode = call(cmd, stdin=stdin,
                                           environment=environment)
        except OSError as error:
            if 'water' in str(error):
                error = 'Water aligner is not installed or not configured properly'
            elif 'exonerate' in str(error):
                error = 'Exonerate aligner is not installed or not configured properly'
            elif 'blast+' in str(error):
                error = 'Blast+ aligner is not installed or not configured properly'
            elif 'lucy' in str(error):
                error = 'Lucy sequence cleaner is not installed or not configured properly'
            elif 'mdust' in str(error):
                error = 'Mdust is not installed or not configured properly'
            elif 'trimpoly' in str(error):
                error = 'Trimpoly sequence cleaner is not installed or not configured properly'
            elif 'sputnik' in str(error):
                error = 'Sputnik microsatellite searcher is not installed or not configured properly'
            elif 'estcan' in str(error):
                error = 'Estcan is not installed or not configured properly'
            raise OSError(error)


        for key, value in runner_data['input'].items():
            for key2, value2 in value.items():
                if key2 == 'fhands':
                    for fhand in value2:
                        fhand.close()

        # there is a error
        if retcode:
            ignore_error = False
            if 'ignore_stderrs' in runner_data:
                for error in runner_data['ignore_stderrs']:
                    if error in stderr:
                        ignore_error = True
            if ignore_error:
                try:
                    print_name = sequence.name
                except AttributeError:
                    print_name = ''

                logging.warning(print_name + ':' + stderr)
            else:
                raise RuntimeError('Problem running ' + tool + ': ' + stdout +
                                   stderr)

        # Now we are going to make this list with the files we are going to
        # return
        returns = {}
        for key, values in runner_data['output'].items():
            #print key, fhand
            if values['option'] == STDOUT:
                fhands = StringIO.StringIO(stdout)
            else:
                fhands = [DisposableFile(fpath) for fpath in values['fpaths']]
                if len(fhands) == 1:
                    fhands = fhands[0]

            returns[key] = fhands
        return returns
    return run_cmd_for_sequence

def _which_binary(binary):
    'It return the full path of the binary if exists'
    stdout = subprocess.PIPE
    process = subprocess.Popen(['which', binary], stdout=stdout)
    stdout = process.communicate()[0]

    if stdout and stdout[0] == '/':
        return stdout.strip()
    else:
        return None

_EXTERNAL_BIN_DIR = None
def call(cmd, environment=None, stdin=None, raise_on_error=False,
         stdout=None, stderr=None, log=False, add_ext_dir=True):
    'It calls a command and it returns stdout, stderr and retcode'
    def subprocess_setup():
        ''' Python installs a SIGPIPE handler by default. This is usually not
        what non-Python subprocesses expect.  Taken from this url:
        http://www.chiark.greenend.org.uk/ucgi/~cjwatson/blosxom/2009/07/02#
        2009-07-02-python-sigpipe'''
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    binary_name = cmd[0]

    if add_ext_dir:
        global _EXTERNAL_BIN_DIR
        if not _EXTERNAL_BIN_DIR:
            ext_dir = get_franklin_ext_dir()
            arch    = platform.architecture()[0]
            system  = platform.system().lower()
            _EXTERNAL_BIN_DIR = os.path.join(ext_dir, 'bin', system, arch)
        binary = os.path.join(_EXTERNAL_BIN_DIR, binary_name)
        cmd[0] = binary

    # emboss binaries need acd files and its environment variable to find them
    if binary_name in ('water', 'est2genome'):
        acd_dir = os.path.join(_EXTERNAL_BIN_DIR, 'emboss_data')
        if environment is None:
            environment = {}
        environment['EMBOSS_ACDROOT'] =  acd_dir
        environment['EMBOSS_DATA']    =  acd_dir


    if stdin is None:
        pstdin = None
    else:
        pstdin = subprocess.PIPE
    if stdout is None:
        stdout = subprocess.PIPE
    if stderr is None:
        stderr = subprocess.PIPE
    #we want to inherit the environment, and modify it
    if environment is not None:
        new_env = {}
        for key, value in os.environ.items():
            new_env[key] = value
        for key, value in environment.items():
            new_env[key] = value
        environment = new_env

    if log:
        logger = logging.getLogger('franklin')
        logger.info('Running command: ' + ' '.join(cmd))

    try:
        process = subprocess.Popen(cmd, stdout=stdout, stderr=stderr,
                                   env=environment, stdin=pstdin,
                                   preexec_fn=subprocess_setup)
    except OSError:
        #if it fails let's be sure that the binary is not on the system
        binary = _which_binary(binary_name)
        if binary is None:
            if 'makeblastdb' in binary_name:
                raise OSError('This program is not installed or not configured properly: ' + cmd[0] + ', part of Blast')
            else:
                raise OSError('This program is not installed or not configured properly: ' + cmd[0])
        #let's try with an absolute path, sometimes works
        cmd.pop(0)
        cmd.insert(0, binary)

        process = subprocess.Popen(cmd, stdout=stdout, stderr=stderr,
                                       env=environment, stdin=pstdin,
                                       preexec_fn=subprocess_setup)

    if stdin is None:
        stdout_str, stderr_str = process.communicate()
    else:
        stdout_str, stderr_str = process.communicate(stdin)
    retcode = process.returncode

    if stdout != subprocess.PIPE:
        stdout.flush()
        stdout.seek(0)
    if stderr != subprocess.PIPE:
        stderr.flush()
        stderr.seek(0)

    if raise_on_error and retcode:
        if stdout != subprocess.PIPE:
            stdout_str = open(stdout.name).read()
        if stderr != subprocess.PIPE:
            stderr_str = open(stderr.name).read()
        msg = 'Error running command: %s\n stderr: %s\n stdout: %s' % \
                                                (' '.join(cmd), stderr_str,
                                                 stdout_str)
        raise RuntimeError(msg)

    return stdout_str, stderr_str, retcode

def b2gpipe_runner(blast, annot_fpath, dat_fpath=None, prop_fpath=None,
                   java_conf=None):
    'It runs b2gpipe'
    java_dir = guess_jar_dir('blast2go.jar', java_conf)
    b2g_bin = os.path.join(java_dir, 'blast2go.jar')
    tempdir = NamedTemporaryDir()
    out_basename = os.path.join(tempdir.name, 'out')

    cmd = java_cmd(java_conf)
    cmd.extend(['-jar', b2g_bin, '-in', blast.name, '-out', out_basename, '-a'])

    if prop_fpath is None:
        prop_fpath = os.path.join(java_dir, 'b2gPipe.properties')
    cmd.extend(['-prop', prop_fpath])

    if dat_fpath:
        cmd.append('-d')
    logger = logging.getLogger('franklin')
    logger.info('Running blast2go: %s' % ' '.join(cmd))

    call(cmd, raise_on_error=True, add_ext_dir=False)
    shutil.move(out_basename + '.annot', annot_fpath)
    if dat_fpath:
        shutil.move(out_basename + '.dat', dat_fpath)
    tempdir.close()

def java_cmd(java_conf):
    'It returns the java -Xmxim thing'
    cmd = ['java']
    if java_conf is None:
        return cmd

    if 'java_memory' in java_conf and java_conf['java_memory'] is not None:
        cmd.append('-Xmx%im' % int(java_conf['java_memory']))
    return cmd
