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

from biolib.seqvar.seqvariation import (calculate_pic, cap_enzime)


## Filters
def create_high_variable_region_filter(max_variability, window_length):
    'It creates a filter that filters seq_vars from high variable regions'
    def high_variable_region_filter(seq_var):
        'The filter'
        if seq_var is None:
            return False
        (seq_var, context) = seq_var

        #how many snps are in the window?
        seq_var_quantity = len(context)
        seq_var_percent = seq_var_quantity / float(window_length * 2)
        if seq_var_percent > max_variability:
            return True
        else:
            return False
    return high_variable_region_filter

def create_close_to_seqvar_filter(distance):
    '''It returns a filter that filters snps by the proximity to other snps.

    If the seqvar has another seqvar closer than DISTANCE, then this snp is
    filtered out'''
    def close_to_seqvar_filter(seq_var):
        'The filter'
        if seq_var is None:
            return False
        (seq_var, context) = seq_var
        location = seq_var.location
        for seqvar_in_contig in context:
            dist = abs(seqvar_in_contig.location - location)
            if dist == 0 :
                continue
            if dist < distance:
                return False
        return True
    return close_to_seqvar_filter

def create_major_allele_freq_filter(frequency):
    '''It creates a filter in which the seqvar is filtered if the most abundant
    allele frequency is bigger than the given one'''

    def major_allele_freq_filter(seq_var):
        'The filter'
        if seq_var is None:
            return False
        seq_var = seq_var[0]
        alleles     = seq_var.alleles
        read_number = 0
        first_read  = alleles[0]['reads']
        for allele in alleles:
            read_number += allele['reads']
        first_percent = (first_read) / float(read_number)
        if first_percent > frequency:
            return False
        else:
            return True
    return major_allele_freq_filter

def create_pic_filter(min_pic):
    '''This funtion is a factory function that creates a function that look
    for the pic of the seqvar and depending on the pic value it filters the
    seqvar o not'''

    def pic_filter(seq_var):
        'The pic filter'
        if seq_var is None:
            return False
        seq_var = seq_var[0]
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
            return False
        seq_var = seq_var[0]
        if seq_var is None:
            return None
        enzymes = cap_enzime(seq_var, all_enzymes)
        if len(enzymes) != 0:
            return True
        else:
            return False
    return enzymes_filter

def create_allele_number_filter(num_alleles):
    'This function factory creates a filter that filters by allele number'
    def allele_number_filter(seq_var):
        'The filter'
        if seq_var is None:
            return False
        seq_var = seq_var[0]
        if len(seq_var.alleles) < num_alleles:
            return False
        else:
            return True
    return allele_number_filter
def create_seqvar_close_to_limit_filter(max_distance):
    '''This function is a function factory. These functions are filters than
     return true if those seqvars ar not clser to the limita than max_distance
     '''
    def seqvar_close_to_reference_limit(seq_var):
        '''True if the sequence variation is close to the contig limits.

        It checks if the seq_variation is not within the range covered by the
        consensus and if it's close to one of its limits. In both cases it will
        return True.
        '''
        if seq_var is None:
            return None
        seq_var = seq_var[0]
        #where does the sequence variation starts and ends?
        location = seq_var.location

        #where does the consensus starts and ends?
        reference = seq_var.reference
        if reference  is None:
            msg = 'SeqVariation is not associated with an reference sequence'
            raise ValueError(msg)

        reference_len = len(reference)
        if (reference_len - location + 1 > max_distance or
            location  > max_distance):
            return True
        else:
            return False

    return seqvar_close_to_reference_limit

# Cleaners
def create_bad_quality_reads_cleaner(qual_treshold):
    '''This function is a factory function that giving a seqvar removes bad
    quality reads from it'''
    def bad_quality_reads_cleaner(seqvar):
        "The cleaner"
        if seqvar is None:
            return None
        seq_var = seqvar[0]
        alleles = seq_var.alleles
        new_alleles = []
        for allele in alleles:
            new_qual = []
            if 'quality' in allele:
                for qual in allele['quality']:
                    if qual > qual_treshold:
                        new_qual.append(qual)
            # check if allele has changed
            len_new_qual = len(new_qual)
            if len_new_qual == 0:
                continue
            elif len_new_qual != len(allele['quality']):
                new_allele = {}
                new_allele['allele']  = allele['allele']
                new_allele['kind']    = allele['kind']
                new_allele['reads']   = len(new_qual)
                new_allele['quality'] = new_qual
                new_alleles.append(new_allele)
            else:
                new_alleles.append(allele)
        if len(new_alleles) == 0:
            return None
        else:
            return seq_var.copy(alleles=new_alleles)

    return bad_quality_reads_cleaner


