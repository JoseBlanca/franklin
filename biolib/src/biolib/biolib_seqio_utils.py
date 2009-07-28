'''
Created on 2009 uzt 28

@author: peio
'''

import tempfile, StringIO, math
from uuid import uuid4

from biolib.seqs import SeqWithQuality
from biolib.biolib_utils import FileIndex

from Bio import SeqIO

def fasta_str(seq, name, description=None):
    'Given a sequence it returns a string with the fasta'
    fasta_str_ = ['>']
    fasta_str_.append(name)
    if description is not None:
        fasta_str_.append('  ', description)
    fasta_str_.append('\n')
    fasta_str_.append(str(seq).strip())
    fasta_str_.append('\n')
    return ''.join(fasta_str_)

def append_to_fasta(seq, seq_fhand, qual_fhand=None):
    'It appends a SeqWithQuality to a fasta file. fhands must be in w mode'
    name        = seq.name
    sequence    = seq.seq
    description = seq.description
    #write seq fasta
    seq_fasta = fasta_str(sequence, name, description)
    seq_fhand.write(seq_fasta)
    #write qual fasta
    if qual_fhand is not None:
        qual = [str(q) for q in seq.qual]
        qual_fasta = fasta_str(" ".join(qual), name, description)
        qual_fhand.write(qual_fasta)


def _get_seq_name(seq):
    'Given a sequence and its default name it returns its name'
    try:
        name = seq.name
    except AttributeError:
        name = None

    if name is None:
        name  = str(uuid4())
    return name

def write_fasta_file(seqs, fhand_seq, fhand_qual=None):
    '''Given a Seq and its default name it returns a fasta file in a
    temporary file. If the seq is a SeqWithQuality you can ask a qual fasta
    file'''
    try:
        # Is seqs an seq or an iter??
        #pylint:disable-msg=W0104
        seqs.name
        seqs = [seqs]
        #pylint:disable-msg=W0704
    except AttributeError:
        pass

    for seq in seqs:
        name = _get_seq_name(seq)

        if fhand_qual is not None:
            try:
                quality = seq.qual
            except AttributeError:
                msg = 'Sequence must be a SeqWithQuality instance'
                raise AttributeError(msg)
            if quality is not None:
                quality = [str(qual) for qual in quality]
                fhand_qual.write(fasta_str(' '.join(quality), name))
            else:
                raise AttributeError('Quality can not be empty')
        fhand_seq.write(fasta_str(seq, name))

    fhand_seq.flush()
    fhand_seq.seek(0)

    if fhand_qual is not None:
        fhand_qual.flush()
        fhand_qual.seek(0)

def temp_fasta_file(seqs, write_qual=False):
    '''Given a Seq and its default name it returns a fasta file in a
    temporary file. If the seq is a SeqWithQuality you can ask a qual fasta
    file'''
    fhand_seq = tempfile.NamedTemporaryFile(suffix='.fasta')
    if write_qual:
        fhand_qual = tempfile.NamedTemporaryFile(suffix='.fasta')
    else:
        fhand_qual = None
    write_fasta_file(seqs, fhand_seq, fhand_qual)

    if write_qual:
        return fhand_seq, fhand_qual
    else:
        return fhand_seq

def temp_multi_fasta_file(seqs):
    '''It creates a temporari fasta file using a list of sequences'''
    fileh = tempfile.NamedTemporaryFile(suffix='.fasta')
    for seq in seqs:
        name = _get_seq_name(seq)
        fileh.write(fasta_str(seq, name))
    fileh.flush()
    fileh.seek(0)
    return fileh

def create_temp_fasta_files(seq1, seq2):
    'It returns two temporal fasta files.'
    #we create two temp files
    fileh1 = temp_fasta_file(seq1)
    fileh2 = temp_fasta_file(seq2)
    return fileh1, fileh2

def parse_fasta(seq_fhand, qual_fhand=None):
    '''It returns the fasta file content giving a file hanler'''
    seq  = []
    qual = None
    name, description, seq = get_content_from_fasta(seq_fhand)
    if qual_fhand is not None:
        #pylint: disable-msg=W0612
        name_, descrip, qual = get_content_from_fasta(qual_fhand, kind='qual')
    if not seq and not qual:
        return None
    else:
        return SeqWithQuality(seq=seq, qual=qual, name=name,
                          description=description)

def get_content_from_fasta(fhand, kind='seq'):
    '''It returns the seq/qual from a fasta file, it need a fhand'''
    seq         = []
    name        = None
    description = None
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            items = line.split(' ', 1)
            name  = items[0][1:]
            try:
                description = items[1]
            except IndexError:
                description = None
        else:
            seq.append(line)
    if kind == 'seq':
        seq = ''.join(seq)
    else:
        qual = []
        for item in seq:
            qual.extend(item.split())
        qual = [int(qitem) for qitem in qual]
        seq = qual
    return name, description, seq

def guess_seq_file_format(fhand):
    'Given a sequence file it returns its format'
    fhand.seek(0)
    line = fhand.readline()
    if line[0] == '>':
        item = fhand.readline().strip()[0]
        if item.isdigit():
            format_ = 'qual'
        else:
            format_ = 'fasta'
    elif line.split()[0] == 'LOCUS':
        format_ = 'genbank'
    else:
        raise ValueError('Unknown sequence file format for : ' + fhand.name)
    fhand.seek(0)
    return format_

def seqs_in_file(seq_fhand, qual_fhand=None, format=None):
    'It yields a seqrecord for each of the sequences found in the seq file'
    if format is None:
        seq_file_format = guess_seq_file_format(seq_fhand)
    else:
        seq_file_format = format

    seq_iter = SeqIO.parse(seq_fhand, seq_file_format)
    if qual_fhand is None:
        qual_iter = None
    else:
        qual_file_format = guess_seq_file_format(qual_fhand)
        qual_iter = SeqIO.parse(qual_fhand, qual_file_format)

    for seqrec in seq_iter:
        #do we have quality?
        letter_annotations = seqrec.letter_annotations
        qual_name = None
        if qual_iter is not None:
            qual_sec_record = qual_iter.next()
            qual = qual_sec_record.letter_annotations['phred_quality']
            qual_name = qual_sec_record.name
        elif 'phred_quality' in letter_annotations:
            qual = letter_annotations['phred_quality']
        elif 'solexa_quality' in letter_annotations:
            qual = letter_annotations['solexa_quality']
            phred = lambda qual: int(10 * math.log(10**(qual/10.0) + 1, 10))
            qual = [phred(value) for value in qual]
        else:
            qual = None

        seq  = seqrec.seq
        name = seqrec.name

        if qual_name and qual_name != name:
            msg = 'Seqs and quals not in the same order: %s, %s' % (name ,
                                                                    qual_name)
            raise RuntimeError(msg)
        description = seqrec.description
        annotations = seqrec.annotations
        yield SeqWithQuality(seq=seq, qual=qual, name=name,
                            description=description, annotations=annotations)

class FileSequenceIndex(object):
    'It indexes sequence files and it returns seq records'
    def __init__(self, fhand, format=None):
        '''It creates the index.

        If the format is not given it will be guessed
        '''
        if format is None:
            format = guess_seq_file_format(fhand)
        self._format = format
        patterns = {
                    'fasta': {'item_start_patterns':['^>'],
                              'key_patterns':['^>([^ \t\n]+)']},
                    'qual': {'item_start_patterns':['^>'],
                              'key_patterns':['^>([^ \t\n]+)']},
                    }
        item_start_patterns = patterns[format]['item_start_patterns']
        key_patterns        = patterns[format]['key_patterns']
        self._index = FileIndex(fhand, item_start_patterns=item_start_patterns,
                                key_patterns=key_patterns)
    def __getitem__(self, name):
        'It returns a sequence record'
        seq_content = StringIO.StringIO(self._index[name])
        if self._format in ('fasta', 'qual'):
            if self._format == 'fasta':
                name, description, seq = get_content_from_fasta(seq_content)
                return SeqWithQuality(name=name, description=description,
                                      seq=seq)
            else:
                name, description, seq = get_content_from_fasta(seq_content,
                                                                kind='qual')
                return SeqWithQuality(name=name, description=description,
                                      qual=seq)
