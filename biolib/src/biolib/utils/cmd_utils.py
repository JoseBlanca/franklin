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


from biolib.utils.seqio_utils import temp_fasta_file, parse_fasta
from biolib.utils.misc_utils import NamedTemporaryDir

import subprocess, signal, tempfile, os
import StringIO, logging

def call(cmd, environment=None, stdin=None, raise_on_error=False):
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
    #we want to inherit the environment, and modify it
    if environment is not None:
        new_env = {}
        for key, value in os.environ.items():
            new_env[key] = value
        for key, value in environment.items():
            new_env[key] = value
        environment = new_env
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=environment,
                               stdin=pstdin, preexec_fn=subprocess_setup)
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
    'blast': {'parameters': {'database' :{'required':True,  'option': '-d'},
                             'program'  :{'required':True,  'option':'-p'},
                             'expect'   :{'default': 0.0001,'option': '-e'},
                             'nhitsv'   :{'default': 20,    'option':'-v'},
                             'nhitsb'   :{'default': 20,    'option':'-b'},
                             #'megablast':{'default':'T',  'option':'-n'},
                             'alig_format': {'default':7, 'option':'-m'}
                            },
              'output':{'blast':{'option':STDOUT}},
              'input':{'option':'-i'},
              'ignore_stderrs': ['Karlin-Altschul']
              },
    'blast+': {'parameters': {'database' :{'required':True,  'option': '-db'},
                    'program'  :{'required':True,  'option':'-p'},
                    'expect'   :{'default': 0.0001,'option': '-evalue'},
                    'nhitsv'   :{'default': 20,    'option':'-num_descriptions'},
                    'nhitsb'   :{'default': 20,    'option':'-num_alignments'},
                    'alig_format': {'default':5, 'option':'-outfmt'}
                            },
              'output':{'blast+':{'option':'-out'}},
              'input':{'option':'-query'},
              'ignore_stderrs': ['Karlin-Altschul']
              },
    'seqclean_vect':{'parameters':{'vector_db':{'required':True,'option':'-v'},
                                 'no_trim_end':{'default':None, 'option':'-N'},
                                'no_trash_low':{'default':None, 'option':'-M'},
                                   'no_trim_A':{'default':None, 'option':'-A'},
                                     'no_dust':{'default':None, 'option':'-L'},
                                 'min_seq_len':{'required':True}},
                     'output':{'sequence':{'option':'-o', 'files':['seq']}},
                     'input':{'option':ARGUMENT, 'arg_before_params':True}
                    },
    'mdust':{'parameters':{'mask_letter':{'default':'L', 'option' : '-m'},
                          'cut_off'    :{'default':'25', 'option':'-v' }},
             'output':{'sequence':{'option':STDOUT}},
             'input':{'option':ARGUMENT, 'arg_before_params':True}
            },
    'trimpoly':{'parameters':{'min_score':{'option':'-s'},
                              'end':{'option':'-e'},
                              'incremental_dist':{'option':'-l'},
                              'fixed_dist':{'option':'-L'},
                              'only_n_trim':{'option':'-N'},
                              'ntrim_above_percent':{'option':'-n'}
                              },
             'output':{'sequence':{'option':STDOUT}},
             'input':{'option':STDIN}
                },
    'exonerate':{'parameters':{'target':{'required':True, 'option':'--target'},
       'show_vulgar':{'default':'False', 'option':'--showvulgar'},
       'show_alignment':{'default':'False', 'option':'--showalignment'},
     'how_options':{'default':"cigar_like:%S %ql %tl %ps\n", 'option':'--ryo'}},
                 'output':{'exonerate':{'option':STDOUT}},
                 'input' :{'option':'--query'}
                 },
    'lucy':{'parameters':{
                      'cdna'   :{'option':'-c', 'default':None},
                      'keep'   :{'option':'-k', 'default':None},
                      'bracket':{'option':'-b', 'default':[10, 0.02]},
                      'window' :{'option':'-w', 'default':[50, 0.08, 10, 0.3]},
                      'error'  :{'option':'-e', 'default':[0.015, 0.015]},
                      'vector' :{'option':'-vector'}
                      },

            'input':{'option': ARGUMENT,  'arg_before_params':True,
                     'files':['seq', 'qual']},
            'output':{'sequence':{'option': '-output','files':['seq', 'qual']}}
            },
#    'lucy.py':{'parameters':{
#                      'cdna'   :{'option':'-c', 'default':None},
#                      'keep'   :{'option':'-k', 'default':None},
#                      'bracket':{'option':'-b', 'default':[10, 0.02]},
#                      'window' :{'option':'-w', 'default':[50, 0.08, 10, 0.3]},
#                      'error'  :{'option':'-e', 'default':[0.015, 0.015]},
#                      'vector' :{'option':'-vector'}
#                      },
#
#               'input':{'inseq': {'option': '--seq', 'file':'seq'},
#                        'inqual': {'option': '--inqual', 'file':'qual'}},
#               'output':{'sequence':{'option': '--outseq', 'file':'seq'},
#                         'quality' :{'option': '--outqual', 'file':'qual'}},
#            },
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


def create_runner(kind, bin_=None, parameters=None, multiseq=False,
                  environment=None):
    ''''It creates a runner class.

    The runner will be able to run a binary program for different sequences.
    kind is the type of runner (blast, seqclean, etc)
    if multiseq is True the runner will expect list or iterator of sequences.
    '''
    runner_def = RUNNER_DEFINITIONS[kind]
    bin_ = _get_aligner_binary(kind)
    if parameters is None:
        parameters = {}
    if environment is None:
        environment = {}
    cmd_param = _process_parameters(parameters, runner_def['parameters'])

    def run_cmd_for_sequence(sequence):
        'It returns a result for the given sequence or sequences'
        #parameters should be in the scope because some tempfile could be in
        #there. In some pythons this has been a problem.
        assert type(parameters)
        # Do we nedd quality for the input??
        input_ = runner_def['input']

        if multiseq:
            seqs = sequence
        else:
            seqs = [sequence]

        if 'files' in input_ and 'qual' in input_['files']:
            seq_fhand, qual_fhand = temp_fasta_file(seqs=seqs,
                                                    write_qual=True)
            qual_fname = qual_fhand.name
        else:
            seq_fhand = temp_fasta_file(seqs=seqs)

        cmd = [bin_]

        #is the sequence a file-like object?
        seq_fname = seq_fhand.name

        #we add the parameter that informs the program where is the input seq
        if input_['option'] == ARGUMENT:
            if 'files' in input_ and 'qual' in  input_['files']:
                input_cmd = [seq_fname, qual_fname]
            else:
                input_cmd = [seq_fname]

        else: # option
            input_cmd = [input_['option'], seq_fname]

        if input_['option'] == STDIN:
            stdin = seq_fhand.read()
            cmd.extend(cmd_param)
        else:
            if 'arg_before_params' in input_:
                cmd.extend(input_cmd)
                cmd.extend(cmd_param)
            else:
                cmd.extend(cmd_param)
                cmd.extend(input_cmd)
        # check if output is passed as a
            stdin = None

        outputs = runner_def['output']
        # Now we are going to construct the output parameters section
        # We have to check if we have different output options, Output is a list
        # of dict with the different output options
        output_cmd    = []
        #Which are the fhands that will be populated by the external cmd?
        output_fhands = {}
        for key, output in outputs.items():
            option = output['option']
            if 'files' not in output:
                files = ['seq']
            else:
                files = output['files']
            if option != STDOUT:
                output_cmd.append(option)

            for file_ in files:
                if option == STDOUT:
                    output_fhands[key] = STDOUT
                else:
                    fhand = tempfile.NamedTemporaryFile()
                    if len(files) > 1:
                        if key not in output_fhands:
                            output_fhands[key] = []
                        output_fhands[key].append(fhand)
                    else:
                        output_fhands[key] = fhand
                    output_cmd.append(fhand.name)
            cmd.extend(output_cmd)
        #print ' '.join(cmd)
        stdout, stderr, retcode = call(cmd, stdin=stdin, environment=environment)

        # there is a error
        if retcode:
            ignore_error = False
            if 'ignore_stderrs' in runner_def:
                for error in runner_def['ignore_stderrs']:
                    if error in stderr:
                        ignore_error = True
            if ignore_error:
                try:
                    print_name = sequence.name
                except AttributeError:
                    print_name = ''

                logging.warning(print_name + ':' + stderr)
            else:
                raise RuntimeError('Problem running ' + bin_ + ': ' + stdout +
                               stderr)

        # Now we are going to make this list with the files we are going to
        # return
        returns = {}
        for key, fhand in output_fhands.items():
            #print key, fhand
            if fhand == STDOUT:
                fhand = StringIO.StringIO(stdout)

            returns[key] = fhand
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

def _get_aligner_binary(kind):
    'it returns the binary to run the given aligner'
    if kind == 'blast':
        return 'blast2'
    else:
        return kind

class SeqcleanRunner(object):
    '''Class to run seqclean '''
    def __init__(self, parameters, ):
        '''Initiator

        We need a working temporal directory to work with seqclean. it poutputs
        a lot of files but we only need two to proceed'''
        self._temp_dir = NamedTemporaryDir()
        self._work_dir = self._temp_dir.name()
        os.chdir(self._work_dir)
        self._parameters = parameters

    def process_sequence(self, sequence):
        ''' Here we process the file '''
        fhand_new_seq = tempfile.NamedTemporaryFile()
        fhand_seq_log = tempfile.NamedTemporaryFile()
        fastah = temp_fasta_file(sequence)
        fastah.flush()

        cmd = ['seqclean', fastah.name]
        cmd.extend(self._parameters)
        cmd.extend(['-o', fhand_new_seq.name])
        cmd.extend(['-r', fhand_seq_log.name])

        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('seqclean run time error:', stderr)
        self._temp_dir.close()
        fhand_new_seq.flush()
        fhand_seq_log.flush()
        seq, name, description = parse_fasta(fhand_new_seq)
        return seq

## This two functions have to be ported to the new runner schema
def translate(seq):
    '''It translates the dna sequence to protein. It uses emboss binary
    transeq'''

    translation_binary = 'transeq'

    fasta_fileh = temp_fasta_file(seq)
    cmd = [translation_binary, fasta_fileh.name, '-stdout', '-auto']
    stdout, stderr, retcode = call(cmd)
    if retcode != 0:
        raise RuntimeError(stderr)
    return parse_fasta(stdout)
def get_best_orf(seq, matrix_path=None):
    '''It returns a new seq with the orf '''

    if matrix_path is None:
        raise ValueError('ESTscan need a matrix to be able to work')
    elif not os.path.exists(matrix_path):
        raise OSError('Matrix file not found: ' + matrix_path)

    estscan_binary = '/usr/local/bin/ESTScan'
    fasta_fileh = temp_fasta_file(seq)
    file_orfh = tempfile.NamedTemporaryFile(suffix='.orf')

    cmd = [estscan_binary, '-M', matrix_path, fasta_fileh.name,
           '-t', file_orfh.name]
    stdout, stderr, retcode = call(cmd)

    if retcode :
        raise RuntimeError(stderr)

    stdout    = StringIO.StringIO(stdout)
    orf_dna  = parse_fasta(stdout)[0]
    orf_prot = parse_fasta(file_orfh)[0]
    return orf_dna, orf_prot
