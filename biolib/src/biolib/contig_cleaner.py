'''
Created on 2009 api 27

@author: peio
'''
from biolib.contig import Location, Contig

def contig_strip(contig, strip_number):
    '''It strips the contig's read ends and starts '''
    striped_contig = Contig(consensus = contig.consensus)
    for read in contig:
        read_location = read.location
        mask_start = read.mask.start + strip_number
        mask_end   = read.mask.end   - strip_number
        if mask_start > mask_end :
            mask = (read.mask.start,read.mask.end )
        else:
            mask = (mask_start, mask_end)
        striped_contig.append_to_location(sequence=read.sequence, \
                                          start=read_location.start, \
                                          strand=read_location.strand,\
                                          forward=read_location.forward,
                                          mask=mask, \
                                          masker=read.masker)
    return striped_contig
