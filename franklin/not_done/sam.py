'''
Code not finished. We shoul change the input format to some standar alignment
format


Created on 2011 api 19

@author: peio
'''
from franklin.sam import sam2bam
from itertools import izip
from tempfile import NamedTemporaryFile

def _reads_in_alignment(fhand):
    'It yields the reads of an alignment'
    ref = None
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        type_, seq = line.split('\t')
        if type_ == 'ref':
            ref = seq
        elif type_.startswith('read'):
            if not ref:
                msg  = 'File malformed, needs ref at the begining of the '
                msg += 'aligment'
                raise ValueError(msg)
            yield ref, seq


def _get_cigar(ref, read):
    'it gets the info needed to build the sam alignment: start and cigar'
    cigar_list = []
    for ref_nucl, read_nucl in izip(ref, read):
        if ref_nucl == '-' and read_nucl == '-' or read_nucl in ('*', ' '):
            continue
        elif ref_nucl != '-' and read_nucl == '-':
            cigar_type = 'D'
        elif ref_nucl != '-' and read_nucl != '-':
            cigar_type = 'M'
        elif ref_nucl == '-' and read_nucl != '-':
            cigar_type = 'I'
        else:
            msg = 'Malformed column: ref_nucl:%s and read_nucl:%s' % \
                                                        (ref_nucl, read_nucl)
            raise ValueError(msg)
        cigar_list.append(cigar_type)


    cigar    = []
    old_type = None
    counter  = 0
    for column_type in cigar_list:
        if old_type is None:
            pass
        elif old_type != column_type:
            cigar.append('%d%s' % (counter, old_type))
            counter = 0
        counter += 1
        old_type = column_type
    else:
        cigar.append('%d%s' % (counter, old_type))
    return ''.join(cigar)

def _get_alignment_start(read):
    'it counts where start the read alignment'
    for index, nucl in enumerate(read):
        index += 1
        if nucl != ' ':
            return str(index)

def sam_creator(fhand, out_bam_path, out_ref_path, read_repeats=None):
    '''it creates a sam using an alignment file. The format of the alignment
    file is:
        ref      aggttttataaaacAAAAaattaagtctacagagcaacta
        sample   aggttttataaaacAAA-aattaagtctacagagcaacta
        read1    aggttttataaaacAA-Aaattaagtctacagagcaacta
        read2    aggttttataaaacA-AAaattaagtctacagagcaacta
        read3    aggttttataaaac-AAAaattaagtctacagagcaacta
    '''
    mapq = '250'

    out_sam = NamedTemporaryFile(suffix='.sam')
    header_done = False

    ref_name = 'ref'
    if read_repeats is None:
        read_repeats = 1
    count = 0
    for ref, read in _reads_in_alignment(fhand):
        ref_seq = ref.replace('-', '').replace('*', '').strip()
        if not header_done:
            out_sam.write('@SQ\tSN:%s\tLN:%d\n' % (ref_name, len(ref_seq)))
            header_done = True
        cigar = _get_cigar(ref, read)
        pos   = _get_alignment_start(read)

        for i in range(read_repeats):
            count += 1
            read_name = 'read%d' % count
            flag  = '0'
            rnext = '*'
            pnext = '0'
            tlen  = '0'
            seq   = read.replace('-', '').replace('*', '').strip()
            qual  = '=' * len(seq)
            sam_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (read_name, flag, ref_name, pos, mapq, cigar, rnext, pnext, tlen, seq, qual)
            out_sam.write(sam_line)


    out_sam.flush()
    sam2bam(out_sam.name, out_bam_path)

    ref_fhand = open(out_ref_path, 'w')
    ref_fhand.write('>ref\n%s' % ref_seq)
    ref_fhand.flush()
