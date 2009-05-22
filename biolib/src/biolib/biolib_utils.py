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
    fasta_str_ = ['>']
    fasta_str_.append(name)
    fasta_str_.append('\n')
    fasta_str_.append(str(seq))
    fasta_str_.append('\n')
    return ''.join(fasta_str_)

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

def remove_from_orf(orf_dna, orf_prot, aminoacid='X'):
    ''' It removes an aminoaacid from dna and protein seq'''
    dna  = []
    prot = []
    pos       = 0
    for letter in orf_prot:
        if letter.upper() != aminoacid:
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
def _remove_atributes_to_tag(tag):
    '''It removees atributes to a xml tag '''
    mod_tag = "".join(tag).split(' ')
    if len(mod_tag) >1:
        return mod_tag[0] + '>'
    else:
        return mod_tag[0]

def _get_xml_header(fhand, tag):
    '''It takes the header of the xml file '''
    fhand.seek(0, 2)
    end_file = fhand.tell()

    fhand.seek(0, 0)
    header      = []
    current_tag = []
    listed_tag = '<' + tag + '>'
    while True:
        if end_file <= fhand.tell():
            raise ValueError('End Of File. Tag Not found')
            
        letter = fhand.read(1)
        if letter == '<':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        else:
            header.append(letter)
        if letter == '>':
            mod_tag = _remove_atributes_to_tag(current_tag)
            if listed_tag == mod_tag:
                return  ''.join(header)
            else:
                header.extend(current_tag)
                current_tag = []
            in_tag = False
def _get_xml_tail(fhand, tag):
    '''It takes the tail of the xml file '''
    in_tag = False
    tail = []
    current_tag  = []
    fhand.seek(-1, 2)
    listed_tag = list('</'+ tag +'>')
    listed_tag.reverse()
    while True:
        if fhand.tell() == 0:
            raise ValueError('Start Of File. Tag Not found')
        letter = fhand.read(1)
        if letter == '>':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        else:
            tail.append(letter)
        
        if letter == '<':
            if current_tag == listed_tag:
                tail.reverse()
                return "".join(tail)
            else:
                tail.extend(current_tag)
                current_tag = []
            in_tag = False 
        
        fhand.seek(-2, 1)
    
def xml_itemize(fhand, tag):  
    '''It takes a xnl file and it chunks it by the given key. It adds header if 
    exists to each of the pieces. It is a generator'''
    fhand.seek(0, 2)
    end_file = fhand.tell()
    
    header = _get_xml_header(fhand, tag)
    tail   = _get_xml_tail(fhand, tag)
    section      = []
    current_tag  = []
    listed_tag_s = '<' + tag + '>'
    listed_tag_e = '</' + tag + '>'
    in_tag, in_section = False, False
    
    fhand.seek(0, 0)
    
    while True:
        if end_file <= fhand.tell():
            break
        letter = fhand.read(1)
        if letter == '<':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        if in_section:
            section.append(letter)            
        if letter == '>':
            if listed_tag_s == _remove_atributes_to_tag(current_tag):
                in_section = True
                section.extend(current_tag)
            elif listed_tag_e == _remove_atributes_to_tag(current_tag):
                yield  header + "".join(section) + tail
                section = []
                in_section = False
            in_tag = False
            current_tag = []
        
    
    
    
 
if __name__ == '__main__':
    pass