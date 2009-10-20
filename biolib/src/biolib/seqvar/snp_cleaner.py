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

from biolib.seqvar.seqvariation import (cap_enzymes, SNP, COMPLEX,
                                        reference_variability,
                                        major_frec_allele_per_library)

#filters
def create_high_variable_region_filter(max_variability, window=None):
    'It creates a filter that filters seq_vars from high variable regions'
    def high_variable_region_filter(snv):
        'The filter'
        if snv is None:
            return False
        (snv, context) = snv
        region_variability = reference_variability(snv, context, window=window)
        if region_variability > max_variability:
            return False
        else:
            return True
    return high_variable_region_filter

def create_close_to_seqvar_filter(distance):
    '''It returns a filter that filters snv by the proximity to other snvs.

    If the snv has another snv closer than DISTANCE, then this snv is
    filtered out'''
    def close_to_snv_filter(snv):
        'The filter'
        if snv is None:
            return False
        (snv, context) = snv
        location = snv.location
        for snv_in_contig in context:
            dist = abs(snv_in_contig.location - location)
            if dist == 0 :
                continue
            if dist < distance:
                return False
        return True
    return close_to_snv_filter

def create_snv_close_to_limit_filter(max_distance):
    '''This function is a function factory. These functions are filters than
     return true if those seqvars ar not clser to the limita than max_distance
     '''
    def snv_close_to_reference_limit(snv):
        '''True if the sequence variation is close to the contig limits.

        It checks if the seq_variation is not within the range covered by the
        consensus and if it's close to one of its limits. In both cases it will
        return True.
        '''
        if snv is None:
            return False
        (snv, context) = snv
        #where does the sequence variation starts and ends?
        location = snv.location

        #where does the consensus starts and ends?
        reference = snv.reference
        if reference  is None:
            msg = 'SeqVariation is not associated with an reference sequence'
            raise ValueError(msg)

        reference_len = len(reference)
        if ((reference_len - location  > max_distance) and
            (location  > max_distance)):
            return True
        else:
            return False

    return snv_close_to_reference_limit


#cleaners
def create_major_allele_freq_filter(frequency, libraries=None):
    '''It creates a cleaner in which the alleles from each library are removed
    if the most abundant allele frequency is bigger than the given one'''
    def major_allele_freq_cleaner(snv):
        'The cleaner'
        if snv is None:
            return False
        snv = snv[0]
        for library_info in snv.per_lib_info:
            if libraries and library_info['library'] not in libraries:
                continue
            alleles = library_info['alleles']
            read_number = 0
            first_read  = alleles[0]['reads']
            for allele in alleles:
                read_number += allele['reads']
            first_percent = (first_read) / float(read_number)
            if first_percent < frequency:
                return True
        return False

    return major_allele_freq_cleaner

def create_bad_quality_reads_cleaner(qual_treshold):
    '''This function is a factory function that giving a seqvar removes bad
    quality reads from it'''
    def bad_quality_reads_cleaner(snv):
        "The cleaner"
        if snv is None:
            return None
        (snv, context) = snv
        new_library_alleles = []
        for library_info in snv.per_lib_info:
            alleles = library_info['alleles']
            new_alleles = []
            for allele in alleles:
                new_qual = []
                if 'qualities' in allele:
                    for qual in allele['qualities']:
                        if qual > qual_treshold or qual is None:
                            new_qual.append(qual)
                else:
                    new_alleles.append(allele)
                    continue

                len_new_qual = len(new_qual)
                if len_new_qual == 0:
                    continue
                elif len_new_qual != len(allele['qualities']):
                    new_allele = {}
                    new_allele['allele']  = allele['allele']
                    new_allele['kind']    = allele['kind']
                    new_allele['reads']   = len(new_qual)
                    new_allele['qualities'] = new_qual
                    new_alleles.append(new_allele)
                else:
                    new_alleles.append(allele)

            # if we have not remove all the alleles, we add the allele to this
            # library
            if len(new_alleles) != 0:
                new_library = {}
                if 'library' in library_info:
                    new_library['library'] = library_info['library']
                new_library['alleles'] = new_alleles
                new_library_alleles.append(new_library)

        if new_library_alleles:
            return (snv.copy(per_lib_info=new_library_alleles), context)

    return bad_quality_reads_cleaner

def create_allele_number_filter(num_alleles, ignore_diffs_with_ref=True,
                                libraries=None):
    '''This function factory creates a mapper that cleans librarys of the snv
    that have less than num alleles'''
    def allele_number_cleaner(snv):
        'The cleaner'
        if snv is None:
            return False
        snv = snv[0]
        for library_info in snv.per_lib_info:
            if libraries and library_info['library'] not in libraries:
                continue
            alleles = library_info['alleles']
            if ((not ignore_diffs_with_ref and len(alleles) == 1 and
                 alleles[0]['kind'] == SNP) or
                (len(alleles) >= num_alleles)):
                return True

        return False
    return allele_number_cleaner

def create_kind_filter(kinds):
    'This filter filters by kind. It only passes snps or complex snv kinds'
    def kind_filter(snv):
        'The filter'
        if snv is None:
            return False
        snv = snv[0]
        if snv.kind in kinds:
            return True
        return False
    return kind_filter

def create_read_number_cleaner(num_reads):
    '''This function factory remove alleles that have more than num_reads
    reads'''
    def read_number_cleaner(snv):
        'The cleaner'
        if snv is None:
            return None
        (snv, context) = snv
        new_library_alleles = []
        for library_info in snv.per_lib_info:
            alleles = library_info['alleles']
            new_alleles = []
            for allele in alleles:
                if allele['reads'] > num_reads:
                    new_alleles.append(allele)

            # if we have not remove all the alleles, we add the allele to this
            # library
            if len(new_alleles) != 0:
                new_library = {}
                if 'library' in library_info:
                    new_library['library'] = library_info['library']
                new_library['alleles'] = new_alleles
                new_library_alleles.append(new_library)
        if new_library_alleles:
            return (snv.copy(per_lib_info=new_library_alleles), context)
    return read_number_cleaner

def create_alleles_n_cleaner():
    '''This function factory removes alleles that composed by N'''
    def read_number_cleaner(snv):
        'The cleaner'
        if snv is None:
            return None
        (snv, context) = snv
        new_library_alleles = []
        for library_info in snv.per_lib_info:
            alleles = library_info['alleles']
            new_alleles = []
            for allele in alleles:
                if allele['allele'].lower() != 'n'*len(allele['allele']):
                    new_alleles.append(allele)

            # if we have not remove all the alleles, we add the allele to this
            # library
            if len(new_alleles) != 0:
                new_library = {}
                if 'library' in library_info:
                    new_library['library'] = library_info['library']
                new_library['alleles'] = new_alleles
                new_library_alleles.append(new_library)
        if new_library_alleles:
            return (snv.copy(per_lib_info=new_library_alleles), context)
    return read_number_cleaner

def create_cap_enzyme_filter(all_enzymes):
    '''This funtion is a factory function that creates a function that look
    if the seqvar is differently afected by some enzymes'''
    def enzymes_filter(snv):
        'The real filter'
        if snv is None:
            return False
        snv = snv[0]
        if snv is None:
            return None
        enzymes = cap_enzymes(snv, all_enzymes)
        if len(enzymes) != 0:
            return True
        else:
            return False
    return enzymes_filter

def create_major_allele_in_lib_frec_filter(frequency):
    '''This function is a factory function that creates a function taht filters
    by the more spreaded allele frequency'''
    def filter_(snv):
        'the filter'
        if snv is None:
            return False
        snv = snv[0]
        if frequency < major_frec_allele_per_library(snv):
            return True
        else:
            return False
    return filter_


