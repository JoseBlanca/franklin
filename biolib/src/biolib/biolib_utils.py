'''
Created on 2009 api 30

@author: peio
'''
import subprocess
import tempfile
from uuid import uuid4
import os.path
import StringIO

def call(cmd):
    'It calls a command and it returns stdout, stderr and retcode'
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    retcode = process.returncode
    return stdout, stderr, retcode


def fasta_str(seq, name):
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
    fileh.write(fasta_str(seq, name))
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


def get_start_end(location):
    '''It accepts an int, Location or tuple and it returns the start, end,
    forward and strand.'''
    #int
    if isinstance(location, int):
        start = location
        end = location
    #tuple
    elif isinstance(location, tuple):
        start = location[0]
        end = location[1]
    #location
    else:
        start = location.start
        end = location.end
    return start, end

def parse_fasta(fhand):
    '''It returns the fasta file content giving a file hanler'''
    seq = []
    for line in fhand:
        line = line.strip()
        if line.startswith('>'):
            items = line.split()
            name  = items[0][1:]
            try:
                description = items[1]
            except IndexError:
                description = None
            continue
        seq.append(line)
    return "".join(seq), name, description

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

def remove_from_orf(orf_dna, orf_prot, aa='X'):
    ''' It removes an aminoaacid from dna and protein seq'''
    dna  = []
    prot = []
    pos       = 0
    for letter in orf_prot:
        if letter.upper() != aa:
            prot.append(letter)
            dna.append(orf_dna[pos:pos + 3])
        pos += 3
    return "".join(dna), "".join(prot)
    
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

#def iprscan_run(seq):
 
if __name__ == '__main__':
    pass