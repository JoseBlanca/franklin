'''
Created on 2009 uzt 28

@author: peio
'''

import tempfile

from Bio import SeqIO

from franklin.seq.seqs import SeqWithQuality
from franklin.seq.readers import seqs_in_file
from franklin.seq.writers import write_seqs_in_file

def xappend_to_fasta(seq, seq_fhand, qual_fhand=None):
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
