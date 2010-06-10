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

import math, re, json

from Bio import SeqIO
from Bio.Alphabet import (Alphabet, SingleLetterAlphabet, ProteinAlphabet,
                          DNAAlphabet)
from Bio.SeqFeature import ExactPosition, FeatureLocation
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from franklin.seq.seqs import (SeqWithQuality, Seq, SeqFeature,
                               create_seq_from_struct)
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
    if line[0] == '{':
        format_ = 'json'
    elif line[:4] in ('SeqW', 'SeqR'):
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
    'It yields a seqrecord for each of the sequences found in the seq file.'

    if format is None:
        format = guess_seq_file_format(seq_fhand)
    seqs =_seqs_in_file(seq_fhand, qual_fhand=qual_fhand, file_format=format)

    if sample_size is None:
        return seqs
    try:
        num_seqs = num_seqs_in_file(seq_fhand, format)
    except NotImplementedError:
        num_seqs = None

    return take_sample(seqs, sample_size, num_seqs)

def _seqs_in_file(seq_fhand, qual_fhand=None, file_format=None):
    'It yields a seqrecord for each of the sequences found in the seq file'
    # look if seq_fhand is a list or not
    seq_fhand.seek(0)
    if qual_fhand is not None:
        qual_fhand.seek(0)

    if file_format == 'repr':
        return _seqs_in_file_with_repr(seq_fhand=seq_fhand)
    if file_format == 'json':
        return _seqs_in_file_with_json(seq_fhand=seq_fhand)
    else:
        return _seqs_in_file_with_bio(seq_fhand=seq_fhand,
                                      file_format=file_format,
                                      qual_fhand=qual_fhand)

def _seqs_in_file_with_json(seq_fhand):
    'It yields all the sequences in json format in a file'

    to_seq = lambda string: create_seq_from_struct(json.loads(string))

    buffer_ = ''
    for line in seq_fhand:
        line = line.rstrip()
        if not line:    #seqs are divided by empty lines
            if buffer:
                yield to_seq(buffer_)
                buffer_ = ''
        else:
            buffer_ += line
    else:
        if buffer_: #the last seq
            yield to_seq(buffer_)

def _seqs_in_file_with_repr(seq_fhand):
    'It yields all the sequences in repr format in a file'

    for seq_chunk in _seq_chunks_in_repr(seq_fhand):
        yield _cast_to_class(seq_chunk)

def _seqs_in_file_with_bio(seq_fhand, file_format, qual_fhand=None):
    '''It yields a seqrecord for each of the sequences found in the seq file
    using biopython'''
    seq_fhand.seek(0)
    if qual_fhand is not None:
        qual_fhand.seek(0)
    #if the format is None maybe the file is empty
    if file_format is None and not seq_fhand.readline():
        raise StopIteration
    seq_fhand.seek(0)
    if qual_fhand is None:
        seq_iter = SeqIO.parse(seq_fhand, BIOPYTHON_FORMATS[file_format])
    else:
        seq_iter = SeqIO.QualityIO.PairedFastaQualIterator(seq_fhand,
                                                           qual_fhand)
    for seqrec in seq_iter:
        #do we have quality?
        letter_annotations = seqrec.letter_annotations

        if 'phred_quality' in letter_annotations:
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
        description = " ".join(seqrec.description.split(' ')[1:])
        annotations = seqrec.annotations
        seqrec = SeqWithQuality(seq=seq, qual=qual, name=name,
                                description=description,
                                annotations=annotations)
        yield seqrec

def _seq_chunks_in_repr(seq_fhand):
    'It yields the reprs for the SeqRecords, one by one'
    buffer_ = ''
    for line in seq_fhand:
        if line[:4] in ('SeqW', 'SeqR'):
            if buffer_:
                yield buffer_
                buffer_ = ''
            buffer_ += line
    #the last sequence
    if buffer_:
        yield buffer_

def _get_dict_key_values(string):
    'Given a string like 1:2, 3:4 it returns the keys and values'

    #print 'string', string
    keys, values = [], []
    for item in _items_in_list(string):
        key, value = list(_items_in_list(item, ':'))
        keys.append(key)
        values.append(value)
    #keys should be all strings
    new_keys = []
    for key in keys:
        new_keys.append(_cast_to_class(key))
    keys = new_keys
    #print 'keys->', keys
    #print 'values->', values
    return keys, values

def _items_in_list(string, split_char=','):
    'Given a string with commas it yields the items in it'
    #the ( add 1 to the depth and the ) remove 1
    #we only return the items at depth 0
    depth = 0
    last_pos = 0 #last position returned
    quotes = ''
    for position, char in enumerate(string):
        #are we inside a string?
        if not quotes and char in ('"', "'"):
            quotes = char
            continue
        elif quotes and char == quotes:
            quotes = ''
            continue
        if quotes:
            continue

        if char in ('(', '[', '{'):
            depth += 1
        elif char in (')', ']', '}'):
            depth -= 1
        elif char == split_char and depth == 0:
            yield string[last_pos: position].strip()
            last_pos = position + 1
    #the last item
    if last_pos < len(string):
        yield string[last_pos: len(string)].strip()

def _item_is_karg(string):
    'It returns True if the given string seems a karg'
    try:
        equal_pos = string.index('=')
    except ValueError:
        return False
    try:
        paren_pos = string.index('(')
    except ValueError:
        return True
    if equal_pos < paren_pos:
        return True
    else:
        return False

def _cast_simple_list(list_repr, class_name):
    'It casts int, string or bool lists'
    list_repr = list_repr.strip()[1:-1]

    if set(('(', '[', '{', ')', ']', '}')).intersection(set(list_repr)):
        return None  #this is not a simple list

    list_ = []
    for item in list_repr.split(','):
        if not item:
            continue
        item = _cast_to_class(item.strip())
        list_.append(item)

    if class_name == 'tuple':
        return tuple(list_)
    else:
        return list_

def _cast_to_class(class_repr):
    'It parses an repr and it returns the data structure'

    #print 'repr ->', class_repr
    class_repr = class_repr.strip()

    if class_repr == '[]':
        return []
    elif class_repr == '{}':
        return {}
    elif class_repr[0] in ('"', "'"):   #string
        return class_repr.strip(class_repr[0])
    elif class_repr[:2] in ('u"', "u'"):   #unicode string
        return unicode(class_repr[1:].strip(class_repr[1]))
    elif class_repr[0] in ('[', '('):
        if class_repr[0] == '[':
            class_name = 'list'
        else:
            class_name = 'tuple'
        # accelerator for simple lists
        simple_list = _cast_simple_list(class_repr, class_name)
        if simple_list is not None:
            return simple_list
        class_content = class_repr[1:-1]
    elif class_repr[0] == '{':
        class_name = 'dict'
        class_content = class_repr[1:-1]
    elif '(' in class_repr:   #non native class
        first_paren_pos = class_repr.index('(')
        last_paren_pos = len(class_repr) - 1
        last_paren_pos = last_paren_pos if class_repr[last_paren_pos] == ')' else last_paren_pos - 1
        class_name = class_repr[:first_paren_pos ].strip()
        class_content =  class_repr[first_paren_pos + 1:last_paren_pos]
    else:
        #a simple element
        if class_repr == 'None':
            return None
        elif class_repr.isdigit():
            return int(class_repr)
        elif class_repr.replace('.', '').isdigit():
            return float(class_repr)
        elif class_repr == 'True':
            return True
        elif class_repr == 'False':
            return False
        else:   #it should be a string
            return class_repr.strip(class_repr[0])

    #now we have to split the class content by the ,
    #also each item can be an arg or a karg
    #print 'class_name ->', class_name
    #print 'class_content ->', class_content

    args, kargs = [], []
    if class_name == 'dict':
        keys, values = _get_dict_key_values(class_content)
        #we cast the values
        values = [_cast_to_class(value) for value in values]
        args = zip(keys, values)
        #print 'dict_args ->', args
    else:
        #print 'items ->', list(_items_in_list(class_content))
        in_kargs = False
        for item in _items_in_list(class_content):
            if not in_kargs and _item_is_karg(item):
                in_kargs = True
            if in_kargs:
                kargs.append(item)
            else:
                args.append(item)

        #we transform the kargs into a dict
        kargs = [arg.split('=', 1) for arg in kargs]

        #we fix the kargs and the args
        kargs = dict([(arg[0].strip(), _cast_to_class(arg[1])) for arg in kargs])
        args = [_cast_to_class(arg) for arg in args]

    classes = {'SeqWithQuality': SeqWithQuality,
               'Seq': Seq,
               'Alphabet': Alphabet,
               'SingleLetterAlphabet': SingleLetterAlphabet,
               'ProteinAlphabet': ProteinAlphabet,
               'ExactPosition': ExactPosition,
               'FeatureLocation': FeatureLocation,
               'DNAAlphabet': DNAAlphabet,
               'SeqFeature': SeqFeature,
               'list': list,
               'tuple': tuple,
               'dict': dict}
    #print 'args->', args
    #print 'kargs->', kargs
    if class_name in ('dict', 'list', 'tuple'):
        return classes[class_name](args)
    else:
        return classes[class_name](*args, **kargs)

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
