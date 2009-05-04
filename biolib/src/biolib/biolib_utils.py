'''
Created on 2009 api 30

@author: peio
'''
import subprocess
import tempfile
from uuid import uuid4
import os.path

def call(cmd):
    'It calls a command and it returns stdout, stderr and retcode'
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    retcode = process.returncode
    return stdout, stderr, retcode


def _fasta_str(seq, name):
    'Given a sequence it returns a string with the fasta'
    fasta_str = ['>']
    fasta_str.append(name)
    fasta_str.append('\n')
    fasta_str.append(str(seq))
    fasta_str.append('\n')
    return ''.join(fasta_str)

def _get_seq_name(seq, name):
    'Given a sequence and its default name it returns its name'
    try:
        return seq.name
    except AttributeError:
        return name 

def temp_fasta_file(seq, name=None):
    '''Given a Seq and its default name it returns a fasta file in a
    temporary file'''
    fileh = tempfile.NamedTemporaryFile(suffix='.fasta')
    if name is None:
        name  = _get_seq_name(seq, str(uuid4()))
    fileh.write(_fasta_str(seq, name))
    fileh.flush()
    return fileh

def create_temp_fasta_files(seq1, seq2):
    'It returns two temporal fasta files.'
    #if the seqs have a name we use it, otherwise we create one
    name1 = _get_seq_name(seq1, 'seq1')
    name2 = _get_seq_name(seq2, 'seq2')
    #we create two temp files
    fileh1 = temp_fasta_file(seq1, name1)
    fileh2 = temp_fasta_file(seq2, name2)
    return fileh1, fileh2


if __name__ == '__main__':
    pass