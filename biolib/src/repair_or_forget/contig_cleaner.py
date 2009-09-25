'''
Created on 2009 api 27

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

from biolib.contig import  Contig
from biolib.seqs import SeqWithQuality
from biolib.locatable_sequence import Location
from biolib.pairwise_alignment import water

def create_contig_read_stripper(length_to_strip=None):
    '''This is a function factory that strips the contigs read  extremes with
    the given length '''
    def contig_strip(contig):
        '''It strips the contig's read ends and starts '''
        if contig is None:
            return None

        if length_to_strip is None:
            raise ValueError('Number of nucleotides to strip needed')
        striped_contig = Contig(consensus = contig.consensus)
        for read in contig:
            read_location = read.location
            mask_start = read.mask.start + length_to_strip
            mask_end   = read.mask.end   - length_to_strip
            if mask_start > mask_end :
                mask = (read.mask.start, read.mask.end )
            else:
                mask = (mask_start, mask_end)
            striped_contig.append_to_location(sequence=read.sequence, \
                                              start=read_location.start, \
                                              strand=read_location.strand,\
                                              forward=read_location.forward,
                                              mask=mask, \
                                              masker=read.masker)
        return striped_contig
    return contig_strip

def create_non_matched_region_stripper():
    '''This is a function factory that creates a function that recalculates the
    mask after aligning consensus with the read'''
    def water_alignment_strip(contig):
        ''' It recalculate de mask after aligning consensus with the read and
        made an intersection between them'''
        if contig is None:
            return None

        consensus = contig.consensus
        new_contig = Contig(consensus=consensus)
        for read in contig:
            read_loc     = read.location

            #with str we get the complemented/reverse sequence if is needed
            seq          = str(read).strip()

            # Alignment of the read with the consensus. The read is only
            # the unmasked region
            align_result = water(consensus, SeqWithQuality(seq, name='seq2'))
            seq_start    = align_result['alignment']['seq2']['start']
            seq_end      = align_result['alignment']['seq2']['end']

            # Now we need to translate the aligned read start and end to the
            # unmasked read coordenates.
            alignment_start_read = read.mask.start + seq_start
            alignment_end_read   = read.mask.start + seq_end

            loc = Location(start=alignment_start_read, end=alignment_end_read, \
                           strand=read.mask.strand, forward=read.mask.forward)

            # Now we have to check with read area are inside the mask and are
            #aligned to the consensus. For that we perform an intersection.
            new_mask = read.mask.intersection(loc)

            # if the intersection doesn't fit, this read is not good for us
            # and we don't use it. If it does fit. we  put inside the contig
            # giving a new mask to the read
            if new_mask is not None:
                mask = (new_mask.start, new_mask.end)
                new_contig.append_to_location(sequence=read.sequence, \
                                              start=read_loc.start, \
                                              strand=read_loc.strand,\
                                              forward=read_loc.forward,
                                              mask=mask, \
                                              masker=read.masker)
        return new_contig
    return water_alignment_strip

def create_read_number_contig_filter(min_read_number):
    '''It is a function factory that return a function that is able to return
    true or false depending on the number of reads that it has'''
    def filter_by_read_number(contig):
        '''It return True if the contig has mora readas than min_read_number
        and false if it has not'''
        if len(contig) > min_read_number:
            return True
        else:
            return False
    return filter_by_read_number









