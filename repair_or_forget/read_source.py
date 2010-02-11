'''
Created on 2009 eka 22

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

import re
from franklin.file_parsers import library_parser, read_clone_library_parser 

class ReadSourceFile(object):
    '''This class '''
    def __init__(self, fhand_library_read_clone):
        '''Initiator '''
        self._fhand_library_read_clone = fhand_library_read_clone
    
    def get_library(self, read):
        '''It returns the library name of the given read'''
        self._fhand_library_read_clone.seek(0)
        for row in read_clone_library_parser(self._fhand_library_read_clone):
            if row['read'].strip() == read:
                return row['library'].strip()
    def get_clone(self, read):
        '''It returns the strain name of the given read'''
        self._fhand_library_read_clone.seek(0)
        for row in read_clone_library_parser(self._fhand_library_read_clone):
            if row['read'].strip() == read:
                return row['clone'].strip()


class ReadSourceRegex(object):
    '''This class accepts a list of regex:library_name parameters and looking 
    on the name of the read it looks if it maches in any of the libraries '''
    def __init__(self, regex_list, naming=None):
        '''It initiates the class.'''
        self._regex_list = regex_list
        self._naming     = naming
        
    def get_library(self, read):
        '''It returns the library name giving a read '''
        for regex_tupla in  self._regex_list:
            regex        = regex_tupla[0]
            library_name = regex_tupla[1]
            if re.match(regex, read):
                return library_name
    
    def get_clone(self, read):
        '''It returns a Naming schema name because this class is going to use 
        only with solexa and 454 (for Now) and those  don't have clones.'''
        raise NotImplementedError('This function is not yet implemented')
    
class ReadSources(object):
    ''' This class accepts a list of read_source classes and it checkes in all 
    of the posibilities to retun a value'''
    
    def __init__(self, read_sources_list):
        'Initiator'
        self._read_sources_list = read_sources_list
    def get_library(self, read):
        'It returns a library checking in the list of classes'
        for read_source_list in self._read_sources_list:
            library = read_source_list.get_library(read)
            if library is not None:
                return library
    def get_clone(self, read):
        'It returns a library checking in the list of classes'
        for read_source_list in self._read_sources_list:
            clone = read_source_list.get_clone(read)
            if clone is not None:
                return clone
              
def get_read_strain(library_name, library_file_list):
    '''It returns the strain name of the given read'''
    for fhand_library in library_file_list:
        for library_dict in library_parser(fhand_library):
            if library_dict['name'] == library_name:
                for cvname, cvtermname, value in library_dict['properties']:
                    if cvtermname == 'strain':
                        return value 
