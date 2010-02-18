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

from Bio import SeqIO

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

class SequenceWriter(object):
    'It writes sequences one by one'
    def __init__(self, fhand, file_format, qual_fhand=None):
        'It inits the class'
        open(fhand.name, 'w')
        self._fhand = open(fhand.name, 'a')

        if qual_fhand is not None:
            open(qual_fhand.name, 'w')
            qual_fhand = open(qual_fhand.name, 'a')
        self._qual_fhand = qual_fhand

        self._format = file_format

    def write(self, sequence):
        'It writes one sequence to the given file'
        format_ = self._format
        if format_ == 'repr':
            self._fhand.write(repr(sequence) + '\n')
        else:
            if ('phred_quality' not in sequence.letter_annotations or
                not sequence.letter_annotations['phred_quality']):
                qual = [30] * len(sequence.seq)
                sequence.letter_annotations['phred_quality'] = qual
            SeqIO.write([sequence], self._fhand, BIOPYTHON_FORMATS[format_])
            if self._qual_fhand and format_ == 'fasta':
                SeqIO.write([sequence], self._qual_fhand, 'qual')

