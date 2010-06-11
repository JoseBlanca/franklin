'''
Created on 2009 uzt 28

@author: peio
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

from franklin.seq.seqs import SeqWithQuality
from franklin.seq.readers import seqs_in_file, guess_seq_file_format
from franklin.seq.writers import write_seqs_in_file

def parse_fasta(seq_fhand, qual_fhand=None):
    '''It returns the fasta file content giving a file hanler'''
    qual = None
    name, description, seq = get_content_from_fasta(seq_fhand)
    if qual_fhand is not None:
        #pylint: disable-msg=W0612
        qual = get_content_from_fasta(qual_fhand, kind='qual')[-1]
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

def seqio(in_seq_fhand, out_seq_fhand, out_format,
          in_qual_fhand=None, out_qual_fhand=None, in_format=None):
    'It converts format of the files'
    if not in_format:
        in_format = guess_seq_file_format(in_seq_fhand)
    if (in_qual_fhand is not None or
        out_qual_fhand is not None or
        in_format in ('repr', 'json', 'pickle') or
        out_format in ('repr', 'json', 'pickle')) :
        seqs = seqs_in_file(seq_fhand=in_seq_fhand,
                            qual_fhand=in_qual_fhand,
                            format=in_format)
        write_seqs_in_file(seqs, seq_fhand=out_seq_fhand,
                           qual_fhand=out_qual_fhand,
                           format=out_format)
    else:
        SeqIO.convert(in_seq_fhand, in_format, out_seq_fhand, out_format)
    out_seq_fhand.flush()
    if out_qual_fhand:
        out_qual_fhand.flush()

def cat(infiles, outfile):
    'It concatenates the given files'
    for infile in infiles:
        if infile is None:
            continue
        infile.seek(0)
        for line in infile:
            outfile.write(line)
