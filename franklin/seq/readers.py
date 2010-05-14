'''This module creates the SeqRecords that will be analyzed by the pipelines.

The functions in this module are SeqRecord generators and SeqRecord mappers.

The generators take an input a sequence file (fasta, fastq, repr or whatever
format) and return a SeqRecord generator from it.
The Readers are classes that have a map method capable of taking a SeqRecord
and adding some annotation. An example of this kind of classes would be a
Reader capaple of taking a bam file and adding the SNPs to the SeqRecords.
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

import math, re

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.utils.itertools_ import take_sample
from franklin.utils.misc_utils import get_fhand

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
                     'qual': 'qual', }

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

def num_seqs_in_file(seq_fhand, format=None):
    'It counts seqs in file. '
    seq_fhand = get_fhand(seq_fhand)
    if format is None:
        format = guess_seq_file_format(seq_fhand)

    if format == 'fasta':
        return count_str_in_file(seq_fhand, '^>')
    elif format == 'repr':
        class_name = SeqWithQuality.__class__.__name__.split('.')[-1]
        return count_str_in_file(seq_fhand, "^%s" % class_name)
    elif 'fastq' in format:
        return _num_seqs_in_fastq(seq_fhand)
    else:
        raise NotImplementedError('I can not count this format: %s' % format)

def count_str_in_file(fhand, regex):
    'it counts the string  in file content'
    regex = re.compile(regex)
    pos_at_start = fhand.tell()
    fhand.seek(0)
    counter = 0
    for line in fhand:
        counter += len(re.findall(regex, line))
    fhand.seek(pos_at_start)
    return counter


def _num_seqs_in_fastq(fhand):
    'it counts seqs in a fastq file '
    pos_at_start = fhand.tell()
    fhand.seek(0)
    counter = 0
    for fastq in FastqGeneralIterator(fhand):
        counter += 1
    fhand.seek(pos_at_start)
    return counter

def seqs_in_file(seq_fhand, qual_fhand=None, format=None, sample_size=None):
    '''It yields a seqrecord for each of the sequences found in the seq file.


    '''
    if format is None:
        format = guess_seq_file_format(seq_fhand)
    seqs =_seqs_in_file(seq_fhand, qual_fhand=qual_fhand, format=format)

    if sample_size is None:
        return seqs
    try:
        num_seqs = num_seqs_in_file(seq_fhand, format)
    except NotImplementedError:
        num_seqs = None

    return take_sample(seqs, sample_size, num_seqs)

def _seqs_in_file(seq_fhand, qual_fhand=None, format=None):
    'It yields a seqrecord for each of the sequences found in the seq file'
    # look if seq_fhand is a list or not
    seq_fhand.seek(0)
    if qual_fhand is not None:
        qual_fhand.seek(0)

    if format == 'repr':
        return _seqs_in_file_with_repr(seq_fhand=seq_fhand)
    else:
        return _seqs_in_file_with_bio(seq_fhand=seq_fhand, format=format,
                                      qual_fhand=qual_fhand)

def _seqs_in_file_with_repr(seq_fhand):
    'It yields all the sequences in repr format in a file'
    from Bio.Alphabet import *
    from Bio.SeqFeature import *
    from franklin.seq.seqs import SeqFeature, Seq

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
    seq_fhand.seek(0)
    if qual_fhand is not None:
        qual_fhand.seek(0)
    #if the format is None maybe the file is empty
    if format is None and not seq_fhand.readline():
        raise StopIteration
    seq_fhand.seek(0)

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
            phred = lambda qual: int(10 * math.log(10 ** (qual / 10.0) + 1, 10))
            qual = [phred(value) for value in qual]
        else:
            qual = None

        seq = seqrec.seq
        seq = Seq(str(seq), seq.alphabet)
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

def guess_seq_type(fhand, format, limit):
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

