'''
Created on 2009 mai 29

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

import urllib2, re

class DatabasesInfo(object):
    ''' It takes db information from a file and it parses it '''
    def __init__(self, source=None):
        ''' Initiator '''
        self.dictionary = None
        if source is None:
            fhand = urllib2.urlopen(
                            "ftp://ftp.geneontology.org/pub/go/doc/GO.xrf_abbs")
        else:
            fhand = open(source, 'r')
        self._parse_abbrs_file(fhand)

    def _parse_abbrs_file(self, fhand):
        '''It parses the file and it puts in a  '''
        db_abbrs = {}
        first = True
        db_dict = {}
        for line in fhand:
            if line[0] == '!' or line.isspace():
                continue
            line = line.strip()
            item, value = line.split(':', 1)
            value = value.strip()
            item  = item.strip()
            if item == 'abbreviation':
                abbrs = line.split(':')[1].strip()
                if first:
                    first = False
                else:
                    db_abbrs[abbrs] = db_dict
                db_dict = {}

            else:
                db_dict[item] = value
        # This one is for the last database
        db_abbrs[abbrs] = db_dict
        self.dictionary = db_abbrs

    def example_id(self, abbrs):
        '''It returns the example id for a abbrs'''
        try:
            return self.dictionary[abbrs]['example_id']
        except KeyError:
            return None
    def database_name(self, abbrs):
        '''It returns the database for a abbrs'''
        return self.dictionary[abbrs]['database']
    def url_syntax(self, abbrs, my_id):
        '''It returns the url_sintax for a abbrs'''
        try:
            url_sintax = self.dictionary[abbrs]['url_syntax']
        except KeyError:
            return None
        return re.sub('\[example_id\]', my_id, url_sintax)
