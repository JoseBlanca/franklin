'''
This library is used to manipulate the sam and bam files. You can merge bams and
add some xtra information to sams

Created on 05/01/2010

@author: peio
'''
import os
from tempfile import NamedTemporaryFile

from biolib.utils.cmd_utils import call
from biolib.utils.seqio_utils import seqs_in_file

PICARDPATH = '/usr/local/biology/picard'

def bam2sam(bampath, sampath=None):
    '''It converts between bam and sam. It sampath is not given, it return
    sam content'''
    cmd = ['samtools', 'view', '-h', bampath]
    if sampath:
        cmd.extend(['-o', sampath])
    sam = call(cmd, raise_on_error=True)
    return sam

def sam2bam(sampath, bampath):
    'It converts between bam and sam.'
    cmd = ['samtools', 'view', '-bth', '-o', bampath, sampath]
    call(cmd, raise_on_error=True)

def bamsam_converter(input_, output_):
    'Converts between sam and bam'
    picard_jar = os.path.join(PICARDPATH, 'SamFormatConverter.jar')
    cmd = ['java', '-Xmx2g', '-jar', picard_jar, 'INPUT=' + input_,
           'OUTPUT=' + output_]
    call(cmd, raise_on_error=True)

def add_header_and_tags_to_sam(sam_fhand, new_sam_fhand):
    '''It adds tags to each of the reads of the alignments on the bam.
    It creates a new bam in other to avoid errors'''

    # first we write the headers
    sam_fhand.seek(0)
    for line in sam_fhand:
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            new_sam_fhand.write(line + '\n')

    # we add a new header. We take this information from the name of the file
    prefix = os.path.basename(sam_fhand.name)
    prefix = prefix.split('.')[:-1]

    rgid_ = []
    append_to_header = []
    for item in prefix:
        try:
            key, value = item.split('_', 1)
            if key.upper() not in ['SM', 'LB', 'PT']:
                continue
        except ValueError:
            continue

        append_to_header.append('%s:%s' % (key.upper(), value))
        rgid_.append(value)
    rgid = "+".join(rgid_)

    # a proper sam with RG needs a SM key
    if not _check_tag_in_RG('SM', append_to_header):
        append_to_header.append('SM:%s' % rgid)


    new_sam_fhand.write('@RG\tID:%s\t' % rgid)
    new_sam_fhand.write("\t".join(append_to_header) + " \n")

    # now the alignments
    sam_fhand.seek(0)
    for line in sam_fhand:
        line = line.strip()
        if not line:
            continue
        if not line.startswith('@'):
            line = line + "\t" + 'RG:Z:%s' % rgid
            new_sam_fhand.write(line + '\n')
    new_sam_fhand.flush()

def _check_tag_in_RG(tag, tags_to_append):
    'It check if the given tag is in the tags to append'
    for tag_to_append in tags_to_append:
        key = tag_to_append.split()[0]
        if key.upper() == tag.upper():
            return True
    return False

def merge_sam(infiles, outfile, reference):
    'It merges a list of sam files'

    #first the reference part of the header
    ref_header = []
    for seq in seqs_in_file(reference):
        name   = seq.name
        length = len(seq)
        ref_header.append(['@SQ', 'SN:%s' % name, 'LN:%d' % length])

    #now the read groups
    headers = set()

    for input_ in infiles:
        input_.seek(0)
        for line in input_:
            line = line.strip()
            if line.startswith('@SQ'):
                continue
            elif line.startswith('@'):
                headers.add(tuple(line.split()))
            else:
                break

    #join and writhe both header parts
    ref_header.extend(headers)
    for header in ref_header:
        outfile.write('\t'.join(header))
        outfile.write('\n')

    #the non header parts
    for input_ in infiles:
        input_.seek(0)
        for line in input_:
            if line.startswith('@'):
                continue
            #outfile.write('\t'.join(line.split()))
            #outfile.write('\n')
            outfile.write(line)
    outfile.flush()

def sort_bam_sam(infile, outfile, sort_method='coordinate'):
    'It sorts a bam file using picard'
    picard_sort_jar = os.path.join(PICARDPATH, 'SortSam.jar')
    cmd = ['java', '-Xmx2g', '-jar', picard_sort_jar, 'INPUT=' +  infile,
           'OUTPUT=' + outfile, 'SORT_ORDER=' + sort_method]
    call(cmd, raise_on_error=True)



