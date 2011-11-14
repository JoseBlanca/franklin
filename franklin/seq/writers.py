'''This module represents the reduce part of our pipelines.

All writes share the same interface. The write method takes a SeqRecord and
appends the part that correspond to that sequence (It can be an annotation or a
feature) to the file.
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

from __future__ import division
import tempfile, json
import cPickle as pickle
import re

from Bio import SeqIO

from franklin.seq.seqs import (get_seq_name, fix_seq_struct_for_json,
                               reverse_complement)
from franklin.seq.readers import (BIOPYTHON_FORMATS, guess_seq_file_format,
                                  seqs_in_file)

class OrthologWriter(object):
    'It writes the orthoolog annotation into a file'
    def __init__(self, fhand):
        'It initiates the class'
        self.fhand = fhand
        header = 'Spp\tseqname\torthologs\n'
        self.fhand.write(header)
        self.num_features = 0

    def write(self, sequence):
        'It does the real write of the features'
        name = get_seq_name(sequence)
        for annot_keys, annot_values in sequence.annotations.items():
            if 'ortholog' in annot_keys:
                self.num_features += 1
                spp = annot_keys.split('-')[0]
                orthologs = ','.join(annot_values)
                line_content = '%s\t%s\t%s\n' % (spp, name, orthologs)
                self.fhand.write(line_content)
        self.fhand.flush()

class OrfWriter(object):
    'It writes the orf annotation into a file'
    def __init__(self, fhand, pep_fhand):
        'It initiates the class'
        self.fhand = fhand
        self.pep_fhand = pep_fhand
        self.num_features = 0

    def write(self, sequence):
        'It does the real write of the features'
        for orf in sequence.get_features(kind='orf'):
            self.num_features += 1
            name = get_seq_name(sequence)
            start = int(str(orf.location.start)) + 1
            end = int(str(orf.location.end)) + 1
            strand = orf.qualifiers['strand']
            seq_content = '>%s_orf_seq start=%d end=%d strand=%s\n%s\n' % \
                          (name, start, end, strand, str(orf.qualifiers['dna']))
            pep_content = '>%s_orf_pep start=%d end=%d strand=%s\n%s\n' % \
                          (name, start, end, strand, str(orf.qualifiers['pep']))
            self.fhand.write(seq_content)
            self.pep_fhand.write(pep_content)
            self.fhand.flush()
            self.pep_fhand.flush()

class SsrWriter(object):
    'It writes the microsatellite annotation into a file'

    def __init__(self, fhand):
        'It initiates the class'
        self.fhand = fhand
        header = 'Seqname\tstart\tend\tlength\tscore\tkind\tunit\tnum repeats\n'
        self.fhand.write(header)
        self.num_features = 0

    def write(self, sequence):
        'It does the real write of the features'
        seq_name = get_seq_name(sequence)
        for ssr in sequence.get_features(kind='microsatellite'):
            self.num_features += 1
            start = int(str(ssr.location.start)) + 1
            end = int(str(ssr.location.end)) + 1
            score = int(ssr.qualifiers['score'])
            kind = ssr.qualifiers['type']
            unit = ssr.qualifiers['unit']
            length = end - start + 1
            num_repeats = length / len(unit)
            self.fhand.write('%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\n' % (seq_name,
                                                        start, end, length,
                                                        score, kind, unit,
                                                        num_repeats))
            self.fhand.flush()

class SamWriter(object):
    'it writes sam file'
    def __init__(self, reference_fhand, reads_fhand, output_fhand,
                 keep_unmapped=True):
        "the initiator"
        self._reference_fhand = reference_fhand
        self._output_fhand    = output_fhand
        self._write_header()
        self._keep_unmapped = keep_unmapped
        format_ = guess_seq_file_format(reads_fhand)
        self._read_index = SeqIO.index(reads_fhand.name, format=format_)

    def _write_header(self):
        'It writes the header of the sam file'
        for seq in seqs_in_file(self._reference_fhand):
            self._output_fhand.write("@SQ\tSN:%s\tLN:%d\n" % (seq.name, len(seq)))
        self._output_fhand.flush()

    def write(self, alignment):
        'It writes each alignment of the sam'
        mapped = alignment['mapped']
        if not mapped and not self._keep_unmapped:
            return

        strand    = '+' if 'strand' not in alignment else alignment['strand']
        read_name = alignment['query_name']
        seqrecord = self._read_index[read_name]
        cigar     = self._clipp_cigar(alignment['cigar'], len(seqrecord))

        if strand == '-':
            seqrecord = reverse_complement(seqrecord)
            cigar = self._reverse_cigar(cigar)

        # check cigar length
        assert len(seqrecord) == self._sum_cigar(cigar)

        flag = self._guess_flag(alignment['mapped'], strand)
        seq       = str(seqrecord.seq)
        try:
            qual  = seqrecord.letter_annotations["phred_quality"]
            qual = [chr(num + 33) for num in qual]
            qual = ''.join(qual)
        except KeyError:
            qual = '*'

        ref_name = alignment['reference_name']
        ref_pos  = str(alignment['position'])
        mapq     = str(alignment['mapq']) if 'mapq' in alignment else '255'

        rnext    = alignment['rnext'] if 'rnext' in alignment else '*'
        pnext    = alignment['pnext'] if 'pnext' in alignment else '0'
        tlen     = alignment['tlen'] if 'tlen' in alignment else '0'
        sam_attrs = [read_name, flag, ref_name, ref_pos, mapq, cigar, rnext,
                     pnext, tlen, seq, qual]

        if 'opt_fields' in alignment:
            opt_fields = '/t'.join(alignment['opt_fields'])
            sam_attrs.append(opt_fields)

        self._output_fhand.write('\t'.join(sam_attrs) + '\n')
        self._output_fhand.flush()


    @staticmethod
    def _guess_flag(mapped, strand):
        'It guess the mapping flag'
        if not mapped:
            return '4'
        elif mapped and strand == '+':
            return '0'
        else:
            return '16'

    @staticmethod
    def _reverse_cigar(cigar):
        'It return the resverse complement of a cigar alignment'
        cigar_items = re.findall('[0-9]+[^0-9]*', cigar)
        rever_cigar = reversed(cigar_items)
        return ''.join(rever_cigar)

    @staticmethod
    def _sum_cigar(cigar):
        'It sums the length of a cigar aligment'
        numbers = []
        for align in re.findall('[0-9]+[^0-9]*', cigar):
            type_ = align[-1]
            num   = int(align[:-1])
            if type_ in ('M', 'I', 'S'):
                numbers.append(num)
        return sum(numbers)

    def _clipp_cigar(self, cigar, seq_length):
        'It adds the soft clipping information if needed to cigar'
        cigar_sum = self._sum_cigar(cigar)
        diff = seq_length - cigar_sum
        if diff != 0:
            cigar += '%dS' % diff
        return cigar


class SequenceWriter(object):
    'It writes sequences one by one'
    def __init__(self, fhand, file_format, qual_fhand=None):
        'It inits the class'
        self.fhand = fhand
        self.qual_fhand = qual_fhand
        self._format = file_format
        self.num_features = 0

    def write(self, sequence):
        'It writes one sequence to the given file'
        self.num_features += 1
        if sequence is None:
            return
        format_ = self._format
        if format_ == 'repr':
            self.fhand.write(repr(sequence) + '\n')
        elif format_ == 'json':
            struct = json.dumps(fix_seq_struct_for_json(sequence.struct,
                                                        alleles_to_string=True))
            self.fhand.write(struct + '\n\n')
        elif format_ == 'pickle':
            string = pickle.dumps(sequence)
            self.fhand.write(string + '\n\n')
        else:
            SeqIO.write([sequence], self.fhand, BIOPYTHON_FORMATS[format_])
            if self.qual_fhand and format_ == 'fasta':
                SeqIO.write([sequence], self.qual_fhand, 'qual')
        self.fhand.flush()
        if self.qual_fhand:
            self.qual_fhand.flush()

def write_seqs_in_file(seqs, seq_fhand, qual_fhand=None, format='fasta',
                       default_quality=25):
    '''It writes the given sequences in the given files.

    The seqs can be an iterartor or a list of Biopython SeqRecords or
    SeqWithQualities'''
    if format == 'fasta':
        _write_fasta_file(seqs, seq_fhand, default_quality,
                          fhand_qual=qual_fhand)
    else:
        writer = SequenceWriter(fhand=seq_fhand, qual_fhand=qual_fhand,
                                file_format=format)
        for seq in seqs:
            writer.write(seq)

def _write_fasta_file(seqs, fhand_seq, default_quality=None, fhand_qual=None):
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
        if seq is None:
            continue
        name = get_seq_name(seq)
        if fhand_qual is not None:
            try:
                quality = seq.letter_annotations["phred_quality"]
            except (AttributeError, KeyError):
                if default_quality:
                    quality = [default_quality] * len(seq)
                else:
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
        fasta_seq = fasta_str(seq, name, desc)
        fhand_seq.write(fasta_seq)

    fhand_seq.flush()
    if fhand_qual is not None:
        fhand_qual.flush()

def fasta_str(seq, name, description=None):
    'Given a sequence it returns a string with the fasta'
    fasta_str_ = ['>']
    fasta_str_.append(name)
    if description:
        fasta_str_.append(' %s' % description)
    fasta_str_.append('\n')
    try:
        fasta_str_.append(str(seq.seq).strip())
    except AttributeError:
        fasta_str_.append(str(seq).strip())
    fasta_str_.append('\n')
    return ''.join(fasta_str_)

def temp_qual_file(seqs):
    'Given a qual seq it return a temporary qual fasta file'
    fhand_qual = tempfile.NamedTemporaryFile(suffix='.qual')
    for seq in seqs:
        if seq is None:
            continue
        name = get_seq_name(seq)
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
    _write_fasta_file(seqs, fhand_seq, fhand_qual=fhand_qual)

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

def create_temp_seq_file(seqs, format):
    'It returns a temporary file with the sequence in it'
    seq_fhand = tempfile.NamedTemporaryFile(suffix='.seq')
    if format == 'qual':
        qual_fhand = tempfile.NamedTemporaryFile(suffix='.qual')
        format = 'fasta'
    else:
        qual_fhand = None
    write_seqs_in_file(seqs, seq_fhand, qual_fhand=qual_fhand, format=format)
    return seq_fhand, qual_fhand
