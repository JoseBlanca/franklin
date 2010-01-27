'''
Command utilities for biolib

This module provides utilities to run external commands into biolib
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


from biolib.utils.seqio_utils import (temp_fasta_file, parse_fasta,
                                      temp_qual_file)
from biolib.utils.misc_utils import NamedTemporaryDir

import subprocess, signal, tempfile, os, itertools
import StringIO, logging, copy

def call(cmd, environment=None, stdin=None, raise_on_error=False,
         stdout=None, stderr=None):
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
    try:
        process = subprocess.Popen(cmd, stdout=stdout, stderr=stderr,
                                   env=environment, stdin=pstdin,
                                   preexec_fn=subprocess_setup)
    except OSError:
        raise OSError('No such file or directory, executable was ' + cmd[0])
    if stdin is None:
        stdout, stderr = process.communicate()
    else:
#        a = stdin.read()
#        print a
#        stdout, stderr = subprocess.Popen.stdin = stdin
#        print stdin.read()
        stdout, stderr = process.communicate(stdin)
    retcode = process.returncode
    if raise_on_error:
        if retcode:
            raise RuntimeError(stderr)
    return stdout, stderr, retcode


# Runner definitions, Define here the parameters of the prgrams you want to
# use with this class
STDOUT   = 'stdout'
ARGUMENT = 'argument'
STDIN    = 'stdin'
RUNNER_DEFINITIONS = {
    'blast': {'binary':'blast2',
              'parameters': {'database' :{'required':True,  'option': '-d'},
                             'program'  :{'required':True,  'option':'-p'},
                             'expect'   :{'default': 0.0001,'option': '-e'},
                             'nhitsv'   :{'default': 20,    'option':'-v'},
                             'nhitsb'   :{'default': 20,    'option':'-b'},
                             #'megablast':{'default':'T',  'option':'-n'},
                             'alig_format': {'default':7, 'option':'-m'}
                            },
              'output':{'blast':{'option':STDOUT}},
              'input':{'sequence':{'option':'-i', 'files_format':['fasta']}},
              'ignore_stderrs': ['Karlin-Altschul']
              },
    'blast+': {'binary':'blast+',
            'parameters': {'database' :{'required':True,  'option': '-db'},
                   'program'  :{'required':True,  'option':'-p'},
                   'expect'   :{'default': 0.0001,'option': '-evalue'},
                   'nhitsv'   :{'default': 20,    'option':'-num_descriptions'},
                   'nhitsb'   :{'default': 20,    'option':'-num_alignments'},
                   'alig_format': {'default':5, 'option':'-outfmt'}
                            },
            'output':{'blast+':{'option':STDOUT}},
            'input':{'sequence':{'option':'-query', 'files_format':['fasta']}},
            'ignore_stderrs': ['Karlin-Altschul']
              },
    'seqclean_vect':{'binary':'seqclean_vect',
                     'parameters':{'vector_db':{'required':True,'option':'-v'},
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
    'estscan':{'binary':'ESTScan',
               'parameters':{ 'matrix':{'option':'-M'},
                             'translate':{'option':'-t', 'default':''}},
               'input':{'sequence':{'option':ARGUMENT,
                                    'arg_before_params':False,
                                    'files_format':['fasta']}},
               'output':{'dna':{'option':STDOUT, 'files_format':['fasta']},
                         'protein':{'option': '-t', 'files_format':['fasta']}}
               }
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
                bin_.append( _param_to_str(value))

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
        fpaths.append(tempfile.mkstemp()[1])
    for fpath in fpaths:
        os.remove(fpath)
    return fpaths

def _prepare_output_files(outputs):
    'It prepares inputs taking into account the format'
    for key, value in outputs.items():
        if 'files_format' in value:
            fhands = _get_mktemp_fpaths(len(value['files_format']))
        else:
            fhands = _get_mktemp_fpaths(1)
        outputs[key]['fpaths'] = fhands

def _build_cmd(cmd_params, runner_def):
    'It bulds the cmd line using  the command definitions'
    inputs = runner_def['input']
    outputs = runner_def['output']
    stdin = None
    bin = runner_def['binary']
    cmd_args_begin = []
    cmd_args_end = []

    for parameters in (inputs, outputs):
        for name, parameter in parameters.items():
            fpaths = parameter['fpaths']
            if parameter['option'] == STDIN:
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

#        print runner_data['input']

        stdout, stderr, retcode = call(cmd, stdin=stdin,
                                       environment=environment)

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
                fhands = [open(fpath) for fpath in values['fpaths']]
                if len(fhands) == 1:
                    fhands = fhands[0]

            returns[key] = fhands
        return returns
    return run_cmd_for_sequence


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
        out_seq_fname =  in_seq_fhand.name


    #we store the output file in a StringIO because the temp dir is going to
    #disapear
    #raw_input()
    result = StringIO.StringIO(open(out_seq_fname).read())

    #we go to the old pwd and we close the temp dir
    os.chdir(old_pwd)
    temp_dir.close()
    return result
