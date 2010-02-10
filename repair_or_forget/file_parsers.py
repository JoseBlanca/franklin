'''
Created on 2009 eka 18

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

import csv

def library_parser(fhand):
    ''' It parses the file and yields a dictionary for each  library'''
    fhand.seek(0)
    library = {}
    for line in fhand:
        line = line.strip()
        if line.startswith('library_definition'):
            #a new library starts
            if library:    #there was a previous library
                yield library
                library = {} #for the following library
            continue
        elif not line:
            continue
        elif line.startswith('format-version'):
            continue
        #if we're here line should be key: value1,value2, value3
        key, values = line.split(':', 1)
        key = key.strip()
        values = [value.strip() for value in values.split(',')]
        #some values should be scalars not lists
        if key == 'name':
            values = values[0]
        if key == 'organism':
            genus, specie = values[0].split(' ')
            library['genus']  = genus.strip()
            library['specie'] = specie.strip()
        elif key == 'type':
            cvname, cvtermname = values[0].split(':')
            library['cvname']     = cvname.strip()
            library['cvtermname'] = cvtermname.strip()
        elif key == 'properties':
            prop_list = []
            for value in values:
                items = value.split(':')
                cvnamep     = items[0].strip()
                cvtermnamep = items[1].strip()
                value       = items[2].strip()
                prop_list.append((cvnamep, cvtermnamep, value))
            library[key] = prop_list
        else:
            library[key] = values
    else:
        #the last library
        yield library
    
def snp_summary_parser(fhand):
    '''It parses the file and yields a dictionary for each  snp'''
    snp = {}
    for line in fhand:
        line = line.strip()
        if line.startswith('snp'):
            
            #a new snp starts
            if snp:    #there was a previous snp
                yield snp
                snp = {} #for the following snp
            continue
        elif not line:
            continue
        elif line.startswith('format-version'):
            continue
        key, values = line.split(':', 1)
        key, values = key.strip(), values.strip()
        snp[key] = values
        if key in ('annotations' , ):
            annot = {}
            values = [value.strip() for value in values.split(',')]
            for value in values:
                prop_type, value_ = value.split(':')
                annot[prop_type.strip()] = value_.strip()
            snp[key] = annot
        
        if key in ('alleles'):
            allele_dic = {}
            values = [value.strip() for value in values.split(';')]
            for value in values:
                allele, reads = value.split(':')
                allele_dic[allele] = [read.strip() for read in reads.split(',')]
            snp[key] = allele_dic
    else:   
        yield snp
        
def read_clone_library_parser(fhand):
    '''It parses clone/read/library file it returns a generator with each 
    line'''
    return  csv.DictReader(fhand, delimiter=',')
        
        
     
    
