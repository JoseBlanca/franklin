'''
This library is used to manipulate the sam and bam files. You can merge bams and
add some xtra information to sams

Created on 05/01/2010

@author: peio
'''
from biolib.utils.cmd_utils import call
from tempfile import NamedTemporaryFile

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

def get_bam_header(bampath):
    'It gets the header of a bam file'
    cmd = ['samtools', 'view', '-H', bampath]
    bamheader = call(cmd, raise_on_error=True)
    return bamheader

def merge_bam(bampath_list, out_bam_path, header=None):
    'It merges two bams adding some lines to the header if needed'

    cmd = ['samtools', 'merge', out_bam_path]
    cmd.extend(bampath_list)

    if header:
        #I need to append the new header to bams header
        samheader = get_bam_header(bampath_list[0])
        samheader += header
        samheader_file = NamedTemporaryFile()
        samheader_file.write(samheader)
        samheader_file.flush()
        # add header option and header file
        cmd.extend(['-h', samheader_file.name])

    call(cmd, raise_on_error=True)

def add_tag_to_bam(bam, tags, newbam):
    '''It adds tags to each of the reads of the alignemnets on the bam.
    It creates a new bam in other to avoid errors'''
    tags_to_add = _make_tag_string(tags)

    sam = bam2sam(bam)
    new_sam = []
    for line in sam.split('\n'):
        line = line.strip()
        if not line:
            continue
        if not line.startswith('@'):
            line = line + "\t" + tags_to_add

        new_sam.append(line)
    tempsam = NamedTemporaryFile()
    tempsam.write("\n".join(new_sam))
    tempsam.flush()
    sam2bam(tempsam.name, newbam)


def _make_tag_string(tags):
    '''It converts the tag dictionariin a string that can be appended to each of
    the read on a bam'''
    tags_ = []
    tag_type = {'RG':'Z'}
    for key, value in tags.items():
        tags_.append('%s:%s:%s' % (key, tag_type[key], value))
    return "\t".join(tags_)







