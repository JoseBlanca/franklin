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


from biolib.biolib_utils import temp_fasta_file, NamedTemporaryDir, parse_fasta
import subprocess, signal, tempfile, os
import StringIO

def call(cmd, env=None, stdin=None, ):
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
        pstdin = stdin

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=env, stdin=pstdin,
                               preexec_fn=subprocess_setup)
    if stdin is None:
        stdout, stderr = process.communicate()
    else:
        stdout, stderr = process.communicate(stdin)
    retcode = process.returncode
    return stdout, stderr, retcode


# Runner definitions, Define here the parameters of the prgrams you want to
# use with this class
STDOUT   = True
ARGUMENT = True
RUNNER_DEFINITIONS = {
    'blast': {'parameters': {'database' :{'required':True,  'option': '-d'},
                             'program'  :{'required':True,  'option':'-p'},
                             'expect'   :{'default': 0.0001,'option': '-e'},
                             'nhitsv'   :{'default': 20,    'option':'-v'},
                             'nhitsb'   :{'default': 20,    'option':'-b'},
                             'megablast':{'default':'T',  'option':'-n'},
                             'alig_format': {'default':7, 'option':'-m'}
                            },
              'output':STDOUT,
              'input':{'option':'-i', 'process':temp_fasta_file}
              },
    'seqclean_vect':{'parameters':{'vector_db':{'required':True,'option':'-v'},
                                 'no_trim_end':{'default':None, 'option':'-N'},
                                'no_trash_low':{'default':None, 'option':'-M'},
                                   'no_trim_A':{'default':None, 'option':'-A'},
                                     'no_dust':{'default':None, 'option':'-L'},
                                 'min_seq_len':{'required':True}},
                     'output':'-o',
                     'input':{'option':ARGUMENT,'process':temp_fasta_file,
                            'arg_before_params':True}
                    },
    'mdust':{'parameters':{'mask_letter':{'default':'L', 'option': '-m'}},
             'output':STDOUT,
             'input':{'option':ARGUMENT, 'process':temp_fasta_file,
                      'arg_before_params':True}
            },
    'exonerate':{'parameters':{'target':{'required':True, 'option':'--target'},
         'model' :{'default':'affine:local', 'option': '--model'},
        'show_vulgar':{'default':'False', 'option':'--showvulgar'},
        'show_alignment':{'default':'False', 'option':'--showalignment'},
        'show_options':{'default':"cigar_like:%S %ql %tl\n", 'option':'--ryo'}},
                 'output':STDOUT,
                 'input' :{'option':'--query', 'process':temp_fasta_file}
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

    for param in parameters:
        if param == 'bin':
            continue
        param_opt = parameters_def[param]['option']
        value     = str(parameters[param])
        bin_.extend((param_opt, value))
    return bin_

def create_runner(kind, bin_=None):
    ''''It creates a runner class.

    The runner will be able to run a binary program for different sequences.
    kind is the type of runner (blast, seqclean, etc)
    '''
    runner_def = RUNNER_DEFINITIONS[kind]
    if bin_ is None:
        bin_ = kind
    def init(self, parameters=None):
        'The init method'
        if parameters == None:
            parameters = {}
        self._cmd_param = _process_parameters(parameters,
                                              runner_def['parameters'])

    def get_result(self, sequence):
        'It returns a result for the given sequence'
        #do we have to process the sequence?
        if 'process' in runner_def['input']:
            sequence = runner_def['input']['process'](sequence)
        cmd = [bin_]

        #we add the input parameter to the cmd
        seq_param = sequence

        #is the sequence a file-like object?
        seq_attrs = dir(sequence)
        if 'name' in seq_attrs and 'close' in seq_attrs:
            seq_param = seq_param.name
        #we add the parameter that informs the program where is the input seq
        if runner_def['input']['option'] == ARGUMENT:
            input_cmd = [seq_param]
        else:
            input_cmd = [runner_def['input']['option'], seq_param]
        if 'arg_before_params' in runner_def['input']:
            cmd.extend(input_cmd)
            cmd.extend(self._cmd_param)
        else:
            cmd.extend(self._cmd_param)
            cmd.extend(input_cmd)
        #the output that we have to collect
        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('Problem running ' + bin_ + ': ' + stderr)

        if runner_def['output'] == STDOUT:
            return stdout
    class_name = kind + 'Runner'
    class_dict = {'__init__':init, 'get_result':get_result}
    klass = type(class_name, (object,), class_dict)
    return klass

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
