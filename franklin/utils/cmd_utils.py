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


from franklin.seq.writers import temp_fasta_file, temp_qual_file
from franklin.utils.misc_utils import NamedTemporaryDir

import subprocess, signal, tempfile, os, itertools
import StringIO, logging, copy, shutil

def _locate_file(fpath):
    cmd = ['locate', fpath]
    stdout = call(cmd, raise_on_error=True)[0]
    found_path = None
    for line in stdout.splitlines():
        if fpath in line:
            found_path = line.strip()
            break
    return found_path

def guess_java_install_dir(jar_fpath):
    'It returns the dir path using locate on a jar file'
    java_dir_path = _locate_file(jar_fpath)
    if not java_dir_path:
        msg = '%s was not found in your system and it is required' % jar_fpath
        raise RuntimeError(msg)
    java_dir_path = java_dir_path.replace(jar_fpath, '')
    return java_dir_path

# Runner definitions, Define here the parameters of the programs you want to
# use with this class
STDOUT = 'stdout'
ARGUMENT = 'argument'
STDIN = 'stdin'
RUNNER_DEFINITIONS = {
    'blast': {'binary':'blast2',
              'parameters': {'database' :{'required':True, 'option': '-d'},
                             'program'  :{'required':True, 'option':'-p'},
                             'expect'   :{'default': 0.0001, 'option': '-e'},
                             'nhitsv'   :{'default': 20, 'option':'-v'},
                             'nhitsb'   :{'default': 20, 'option':'-b'},
                             #'megablast':{'default':'T',  'option':'-n'},
                             'alig_format': {'default':7, 'option':'-m'}
                            },
              'output':{'blast':{'option':STDOUT}},
              'input':{'sequence':{'option':'-i', 'files_format':['fasta']}},
              'ignore_stderrs': ['Karlin-Altschul']
              },
#    'blast+': {'binary':'blast+',
#            'parameters': {'database' :{'required':True,  'option': '-db'},
#                   'program'  :{'required':True,  'option':'-p'},
#                   'expect'   :{'default': 0.0001,'option': '-evalue'},
#                  'nhitsv'   :{'default': 20,    'option':'-num_descriptions'},
#                   'nhitsb'   :{'default': 20,    'option':'-num_alignments'},
#                   'alig_format': {'default':5, 'option':'-outfmt'}
#                            },
#            'output':{'blast+':{'option':STDOUT}},
#            'input':{'sequence':{'option':'-query', 'files_format':['fasta']}},
#            'ignore_stderrs': ['Karlin-Altschul']
#              },
    'seqclean_vect':{'binary':'seqclean_vect',
                     'parameters':{'vector_db':{'required':True, 'option':'-v'},
                                 'no_trim_end':{'default':None, 'option':'-N'},
                                'no_trash_low':{'default':None, 'option':'-M'},
                                   'no_trim_A':{'default':None, 'option':'-A'},
                                     'no_dust':{'default':None, 'option':'-L'},
                                 'min_seq_len':{'required':True}},
                     'output':{'sequence':{'option':'-o', 'files':['seq']}},
                     'input':{'sequence':{'option':ARGUMENT,
                                          'arg_before_params':True,
                                          'files_format':['fasta']}}
                    },
    'mdust':{'binary':'mdust',
             'parameters':{'mask_letter':{'default':'L', 'option' : '-m'},
                          'cut_off'    :{'default':'25', 'option':'-v' }},
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
    'exonerate':{'binary':'exonerate',
                 'parameters':{'target':{'required':True, 'option':'--target'},
       'show_vulgar':{'default':'False', 'option':'--showvulgar'},
       'show_alignment':{'default':'False', 'option':'--showalignment'},
     'how_options':{'default':"cigar_like:%S %ql %tl %ps\n", 'option':'--ryo'}},
                 'output':{'exonerate':{'option':STDOUT}},
                 'input' : {'sequence':{'option':'--query',
                                        'files_format':['fasta']}}
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
    'remap':{'binary':'remap',
             'parameters':{'enzymes':{'default':'all', 'option' : '-enzymes'},
                          'sitelen' :{'default':'4', 'option':'-sitelen' },
                          'stdout'  :{'default':'', 'option':'-stdout' },
                          'auto'    :{'default':'', 'option':'-auto' }, },
             'output':{'remap':{'option':'-outfile', 'files':['map']}},
             'input':{'sequence':{'option': '-sequence',
                                  'files_format':['fasta']}}
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

class DisposableFile(file):
    'A file that remove the file when closed'
    def close(self):
        'This close removes the file when called'
        file.close(self)
        os.remove(self.name)

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
    general_cmd_param = _process_parameters(parameters,
                                    RUNNER_DEFINITIONS[tool]['parameters'])

    def run_cmd_for_sequence(sequence):
        'It returns a result for the given sequence or sequences'
        #parameters should be in the scope because some tempfile could be in
        #there. In some pythons this has been a problem.
        runner_data = copy.deepcopy(RUNNER_DEFINITIONS[tool])
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

        stdout, stderr, retcode = call(cmd, stdin=stdin,
                                       environment=environment)


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
    process = subprocess.Popen(['/bin/which', binary], stdout=stdout)
    stdout = process.communicate()[0]

    if stdout and stdout[0] == '/':
        return stdout.strip()
    else:
        return None

def call(cmd, environment=None, stdin=None, raise_on_error=False,
         stdout=None, stderr=None, log=False):
    'It calls a command and it returns stdout, stderr and retcode'
    def subprocess_setup():
        ''' Python installs a SIGPIPE handler by default. This is usually not
        what non-Python subprocesses expect.  Taken from this url:
        http://www.chiark.greenend.org.uk/ucgi/~cjwatson/blosxom/2009/07/02#
        2009-07-02-python-sigpipe'''
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

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
        binary = _which_binary(cmd[0])
        if binary is None:
            raise OSError('The binary was not found: ' + cmd[0])
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
    if raise_on_error and retcode:
        msg = 'Error running command: %s\n stderr: %s\n stdout: %s' % \
                                                (' '.join(cmd), stderr_str,
                                                 stdout_str)
        raise RuntimeError(msg)
    if stdout != subprocess.PIPE:
        stdout.flush()
    if stderr != subprocess.PIPE:
        stderr.flush()
    return stdout_str, stderr_str, retcode

def b2gpipe_runner(blast, annot_fpath, dat_fpath=None, prop_fpath=None,
                   java_memory=None):
    'It runs b2gpipe'
    java_dir = guess_java_install_dir('blast2go.jar')
    b2g_bin = os.path.join(java_dir, 'blast2go.jar')
    tempdir = NamedTemporaryDir()
    out_basename = os.path.join(tempdir.name, 'out')
    cmd = java_cmd(java_memory)
    cmd.extend(['-jar', b2g_bin, '-in', blast.name, '-out', out_basename, '-a'])

    if prop_fpath is None:
        prop_fpath = os.path.join(java_dir, 'b2gPipe.properties')
        cmd.extend(['-prop', prop_fpath])

    if dat_fpath is not None:
        cmd.append('-d')
    call(cmd, raise_on_error=True)
    shutil.move(out_basename + '.annot', annot_fpath)
    if dat_fpath is not None:
        shutil.move(out_basename + '.dat', dat_fpath)
    tempdir.close()

def java_cmd(java_conf):
    'It returns the java -Xmxim thing'
    cmd = ['java']
    if 'java_memory' in java_conf:
        cmd.append('-Xmx%im' % int(java_conf['java_memory']))
    return cmd

def run_repeatmasker_for_sequence(sequence, species='eudicotyledons'):
    '''It returns masked sequence (StrinIO) for the given sequence.
    '''

    #we run repeat masker in a temp dir
    temp_dir = NamedTemporaryDir()
    old_pwd = os.getcwd()
    os.chdir(temp_dir.name)

    #input sequence and output file
    in_seq_fhand = temp_fasta_file(sequence)
    out_seq_fname = in_seq_fhand.name + '.masked'

    #parameters used
    # q         fast search (qq is even faster and less sensitive)
    # species   the species to use (e.g. arabidopsis, eudicotyledons)
    # no_is     do not look for bacterial insertions
    # no_cut    do not excise the repeats found
    # xsmall    repeats lower cased
    cmd = ['RepeatMasker', '-q', '-species', species, '-no_is',
           'no_cut', '-xsmall', in_seq_fhand.name]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('Problem running RepeatMasker: ' + stderr)

    # If there are no repetitive sequences, repreatmasker does not return any
    # file
    if 'No repetitive sequences' in stderr:
        out_seq_fname = in_seq_fhand.name


    #we store the output file in a StringIO because the temp dir is going to
    #disapear
    #raw_input()
    resul_fhand = open(out_seq_fname)
    result = StringIO.StringIO(resul_fhand.read())
    resul_fhand.close()

    #we go to the old pwd and we close the temp dir
    os.chdir(old_pwd)
    temp_dir.close()
    return result
