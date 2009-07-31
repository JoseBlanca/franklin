'''
This module contain cleaners and filters to accept or reject founded seq
variations

Created on 2009 uzt 30

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

from biolib.seqvariation import (CONFIG, SeqVariation, allele_count,
                                 calculate_pic, cap_enzime)
from biolib.contig import slice_to_range


def create_bad_quality_allele_remover(qual_threshold=None,default_quality=None):
    '''This function is a factory function that given a seqvar removes bad
    quality alleles'''


    def remove_bad_quality_alleles(seqvar):
        ''' It removes bad quality alleles given a seqvar'''
        if seqvar is None:
            return None

        if  qual_threshold is None:
            raise ValueError('Quality threshold must be implicit')

        alleles  = seqvar.alleles
        contig   = seqvar.alignment
        location = seqvar.location

        bad_alleles = [] # Alleles we are not going to use
        for allele in alleles:
            # Gaps does not have quality, and it is OK
            if CONFIG.indel_char in allele:
                continue
            quality_allele = _allele_quality(contig, alleles[allele], location,
                                              default_quality)

            # We check it the quality is enought to use
            if quality_allele < qual_threshold:
                bad_alleles.append(allele)
        for allele in bad_alleles:
            del alleles[allele]

        if len(alleles) >1:
            return SeqVariation(alleles=alleles, location=seqvar.location, \
                                alignment=seqvar.alignment)
        else:
            return None
    return remove_bad_quality_alleles

def _allele_quality(contig, reads, location, default_quality):
    ''' It returns the quality of a allele'''
    quality_allele = 0
    for read_num in reads:
        read = contig[read_num]
        # If the location is more than one column we need to obtain de media
        # of all columns involved
        try:
            columns = slice_to_range(location, None)
        except AttributeError:
            start = location.start
            end = location.end
            columns = slice_to_range(slice(start, end + 1), None)


        # We sum all qualities of the allele
        read_allele_quality = _read_allele_quality(columns, read,
                                                   default_quality)
        if read_allele_quality is not None:
            quality_allele += read_allele_quality
    # this is the media of the qualities
    return quality_allele / len(reads)


def _read_allele_quality(columns, read, default_quality):
    '''it returns of the quality of the allele in one read '''
    quality_row = 0
    cols_with_quality = 0
    for column_location in columns:
        #It checks if the read have quality, and if we are giving a
        #default quality for the reads that haven't
        try:
            #The easiest way to get simple columns is using this.
            # Because it complements and reverses if is needed
            nucleotide_quality = read[column_location].qual[0]
        except IndexError:
            #No quality for this read in this position
            continue
        except TypeError:
            if default_quality is None:
                msg = "No Quality in read and no default provided"
                raise ValueError(msg)
            else:
                nucleotide_quality = default_quality
        quality_row += nucleotide_quality
        cols_with_quality += 1

    if cols_with_quality:
        return quality_row / cols_with_quality
    else:
        return None


def create_second_allele_number_filter(number_2allele):
    '''This function is a function factory that creates a function that gives
     True if the second allele of the seqvar is bigger than number2_allele, if
     not is returns false and the sevar is filtered'''

    def second_allele_read_times(seq_var):
        '''It returns True if the second most abundant allele has been read at
        least the given times'''
        if seq_var is None:
            return None
        alleles = seq_var.sorted_alleles()
        if allele_count(alleles[1][1]) >= number_2allele:
            return True
        return False
    return second_allele_read_times

def create_seqvar_close_to_limit_filter(max_distance):
    '''This function is a function factory. These functions are filters than
     return true if those seqvars ar not clser to the limita than max_distance
     '''
    def seqvar_close_to_consensus_limit(seq_variation):
        '''True if the sequence variation is close to the contig limits.

        It checks if the seq_variation is not within the range covered by the
        consensus and if it's close to one of its limits. In both cases it will
        return True.
        '''
        if seq_variation is None:
            return None
        #where does the sequence variation starts and ends?
        seq_var_loc = seq_variation.location
        try:
            #it's a location
            seq_var_start = seq_var_loc.start
            seq_var_end   = seq_var_loc.end
        except AttributeError:
            #it's an int
            seq_var_start = seq_var_loc
            seq_var_end   = seq_var_loc

        #where does the consensus starts and ends?
        contig = seq_variation.alignment
        if contig is None:
            raise ValueError('SeqVariation is not associated with an alignment')
        consensus = contig.consensus
        if consensus is None:
            raise ValueError('The alignment has no consensus')
        #the consensus might be a LocatableSequence or an standard sequence
        try:
            con_loc = consensus.location
            con_start = con_loc.start
            con_end   = con_loc.end
        except AttributeError:
            con_start = 0
            con_end   = len(consensus) - 1

        #Now we can check the limits
        if (seq_var_start < con_start + max_distance or
            seq_var_end   > con_end - max_distance):
            return False
        return True
    return seqvar_close_to_consensus_limit

def create_pic_filter(min_pic):
    '''This funtion is a factory function that creates a function that look
    for the pic of the seqvar and depending on the pic value it filters the
    seqvar o not'''

    def pic_filter(seq_var):
        'The pic filter'
        if seq_var is None:
            return None
        if calculate_pic(seq_var) < min_pic:
            return False
        else:
            return True
    return pic_filter

def create_cap_enzyme_filter(all_enzymes):
    '''This funtion is a factory function that creates a function that look
    if the seqvar is differently afected by some enzymes'''
    def enzymes_filter(seq_var):
        'The real filter'
        if seq_var is None:
            return None
        enzymes = cap_enzime(seq_var, all_enzymes)
        if len(enzymes) != 0:
            return True
        else:
            return False
    return enzymes_filter






