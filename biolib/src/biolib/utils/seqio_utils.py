'''
Created on 2009 uzt 28

@author: peio
'''

import tempfile, StringIO, math
from uuid import uuid4

from biolib.seq.seqs import SeqWithQuality
from biolib.utils.misc_utils import FileIndex

from Bio import SeqIO
from Bio.Seq import UnknownSeq

def fasta_str(seq, name, description=None):
    'Given a sequence it returns a string with the fasta'
    fasta_str_ = ['>']
    fasta_str_.append(name)
    if description is not None:
        fasta_str_.append(' %s' % description)
    fasta_str_.append('\n')
    try:
        fasta_str_.append(str(seq.seq).strip())
    except AttributeError:
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
    #the SeqRecord default
    if name == "<unknown name>":
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
                quality = seq.letter_annotations["phred_quality"]
            except AttributeError:
                msg = 'Sequence must have a phred_quality letter annotation'
                raise AttributeError(msg)
            if quality is not None:
                quality = [str(qual) for qual in quality]
                fhand_qual.write(fasta_str(' '.join(quality), name))
            else:
                raise AttributeError('Quality can not be empty')
        try:
            desc = seq.description
        except AttributeError:
            desc = None
        if desc == "<unknown description>":
            desc = None
        fhand_seq.write(fasta_str(seq, name, desc))

    fhand_seq.flush()
    fhand_seq.seek(0)

    if fhand_qual is not None:
        fhand_qual.flush()
        fhand_qual.seek(0)

def temp_qual_file(seqs):
    'Given a qual seq it return a temporary qual fasta file'
    fhand_qual = tempfile.NamedTemporaryFile(suffix='.qual')
    for seq in seqs:
        name    = _get_seq_name(seq)
        quality = seq.letter_annotations["phred_quality"]

        quality = [str(qual) for qual in quality]
        fhand_qual.write('>%s\n%s\n' % (name , ' '.join(quality)))

    fhand_qual.flush()
    fhand_qual.seek(0)
    return fhand_qual

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

def create_temp_fasta_files(seq1, seq2):
    'It returns two temporal fasta files.'
    #we create two temp files
    fileh1 = temp_fasta_file(seq1)
    fileh2 = temp_fasta_file(seq2)
    return fileh1, fileh2

def parse_fasta(seq_fhand, qual_fhand=None):
    '''It returns the fasta file content giving a file hanler'''
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

def quess_seq_type(fhand, format, limit):
    'It guess seq type: short, long'
    short = 0
    seqs = seqs_in_file(fhand, format=format)
    for i in range(3):
        if len(seqs.next().seq) < limit:
            short += 1
    if short == 3:
        return 'short_seqs'
    else:
        return 'long_seqs'

def guess_seq_file_format(fhand):
    'Given a sequence file it returns its format'
    fhand.seek(0)
    line = fhand.readline()
    if not line:
        return None
    if line[:4] in ('SeqW', 'SeqR'):
        format_ = 'repr'
    elif line[0] == '>':
        item = fhand.readline().strip()[0]
        if item.isdigit():
            format_ = 'qual'
        else:
            format_ = 'fasta'
    elif line[0] == '@' and fhand.name.endswith('.sfastq'):
        format_ = 'fastq'
    elif line.split()[0] == 'LOCUS':
        format_ = 'genbank'
    else:
        raise ValueError('Unknown sequence file format for : ' + fhand.name)
    fhand.seek(0)
    return format_

#the translation between our formats and the biopython formats
BIOPYTHON_FORMATS = {'fasta': 'fasta',
                     'fastq': 'fastq',
                     'sfastq':'fastq',
                     'fastq-sanger': 'fastq',
                     'ifastq': 'fastq-illumina',
                     'fastq-illumina': 'fastq-illumina',
                     'fastq-solexa': 'fastq-solexa',
                     'genbank': 'genbank',
                     'gb': 'genbank',
                     'embl': 'embl',
                     'qual': 'qual',}

def seqs_in_file(seq_fhand, qual_fhand=None, format=None):
    'It yields a seqrecord for each of the sequences found in the seq file'
    if format is None:
        format = guess_seq_file_format(seq_fhand)
    if format == 'repr':
        return _seqs_in_file_with_repr(seq_fhand=seq_fhand)
    else:
        return _seqs_in_file_with_bio(seq_fhand=seq_fhand, format=format,
                                      qual_fhand=qual_fhand)

def _seqs_in_file_with_repr(seq_fhand):
    'It yields all the sequences in repr format in a file'
    from Bio.Alphabet import *
    from Bio.SeqFeature import *
    from Bio.Seq import Seq

    buffer_ = ''
    for line in seq_fhand:
        if line[:4] in ('SeqW', 'SeqR'):
            if buffer_:
                yield eval(buffer_)
                buffer_ = ''
            buffer_ += line
    if buffer_:
        yield eval(buffer_)
    else:
        raise StopIteration

def _seqs_in_file_with_bio(seq_fhand, format, qual_fhand=None):
    '''It yields a seqrecord for each of the sequences found in the seq file
    using biopython'''
    #if the format is None maybe the file is empty
    if format is None and not open(seq_fhand.name).readline():
        raise StopIteration

    seq_iter = SeqIO.parse(seq_fhand, BIOPYTHON_FORMATS[format])
    if qual_fhand is None:
        qual_iter = None
    else:
        qual_iter = SeqIO.parse(qual_fhand, 'qual')

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
        description = " ".join(seqrec.description.split(' ')[1:])
        annotations = seqrec.annotations
        seqrec = SeqWithQuality(seq=seq, qual=qual, name=name,
                                description=description,
                                annotations=annotations)
        yield seqrec

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
                qual = None
            else:
                name, description, qual = get_content_from_fasta(seq_content,
                                                                kind='qual')
                seq = UnknownSeq(length=len(qual))
            seqrec = SeqWithQuality(seq=seq, qual=qual)
            if name is not None:
                seqrec.name = name
            if description is not None:
                seqrec.description = description
            return seqrec

def write_seqs_in_file(seqs, seq_fhand, qual_fhand=None, format='fasta'):
    '''It writes the given sequences in the given files.

    The seqs can be an iterartor or a list of Biopython SeqRecords or
    SeqWithQualities'''
    for seq in seqs:
        if format == 'repr':
            seq_fhand.write(repr(seq) + '\n')
        else:
            if ('phred_quality' not in seq.letter_annotations or
                not seq.letter_annotations['phred_quality']):
                qual = [30] * len(seq.seq)
                seq.letter_annotations['phred_quality'] = qual
            SeqIO.write([seq], seq_fhand, BIOPYTHON_FORMATS[format])
            if qual_fhand and format == 'fasta':
                SeqIO.write([seq], qual_fhand, 'qual')


def seqio(in_seq_fhand, out_seq_fhand, out_format,
          in_qual_fhand=None, out_qual_fhand=None, in_format=None):
    'It converts format of the files'
    if  in_qual_fhand is not None or out_qual_fhand is not None:
        seqs = seqs_in_file(seq_fhand=in_seq_fhand,
                            qual_fhand=in_qual_fhand,
                            format=in_format)
        write_seqs_in_file(seqs, seq_fhand=out_seq_fhand,
                           qual_fhand=out_qual_fhand,
                           format=out_format)
    else:
        SeqIO.convert(in_seq_fhand, in_format, out_seq_fhand, out_format)

def cat(infiles, outfile):
    'It concatenates the given files'
    for infile in infiles:
        if infile is None:
            continue
        infile.seek(0)
        for line in infile:
            outfile.write(line)
