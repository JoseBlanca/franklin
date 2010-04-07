'''
Created on 2009 api 30

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

import tempfile, shutil
import os, re, math
from franklin.seq.seqs import copy_seq_with_quality
import franklin

DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'franklin',
                         'data')

def float_lists_are_equal(list1, list2):
    'Given two lists it checks that all floats are equal'
    for num1, num2 in zip(list1, list2):
        assert floats_are_equal(num1, num2)

def floats_are_equal(num1, num2):
    'Given two numbers it returns True if they are similar'
    if num1 == 0.0:
        if num2 == 0.0:
            return True
        else:
            return False
    log1 = math.log(float(num1))
    log2 = math.log(float(num2))
    return abs(log1 - log2) < 0.01

class NamedTemporaryDir(object):
    '''This class creates temporary directories '''
    def __init__(self):
        '''It initiates the class.'''
        self._name = tempfile.mkdtemp()

    def get_name(self):
        'Returns path to the dict'
        return self._name
    name = property(get_name)
    def close(self):
        '''It removes the temp dir'''
        if os.path.exists(self._name):
            shutil.rmtree(self._name)

    def __del__(self):
        '''It removes de temp dir when instance is removed and the garbaje
        colector decides it'''
        self.close()

def _remove_atributes_to_tag(tag):
    '''It removees atributes to a xml tag '''
    mod_tag = "".join(tag).split(' ')
    if len(mod_tag) >1:
        return mod_tag[0] + '>'
    else:
        return mod_tag[0]

def _get_xml_header(fhand, tag):
    '''It takes the header of the xml file '''
    fhand.seek(0, 2)
    end_file = fhand.tell()

    fhand.seek(0, 0)
    header      = []
    current_tag = []
    listed_tag = '<' + tag + '>'
    while True:
        if end_file <= fhand.tell():
            raise ValueError('End Of File. Tag Not found')

        letter = fhand.read(1)
        if letter == '<':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        else:
            header.append(letter)
        if letter == '>':
            mod_tag = _remove_atributes_to_tag(current_tag)
            if listed_tag == mod_tag:
                return  ''.join(header)
            else:
                header.extend(current_tag)
                current_tag = []
            in_tag = False

def _get_xml_tail(fhand, tag):
    '''It takes the tail of the xml file '''
    in_tag = False
    tail = []
    current_tag  = []
    fhand.seek(-1, 2)
    listed_tag = list('</'+ tag +'>')
    listed_tag.reverse()
    while True:
        if fhand.tell() == 0:
            raise ValueError('Start Of File. Tag Not found')
        letter = fhand.read(1)
        if letter == '>':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        else:
            tail.append(letter)

        if letter == '<':
            if current_tag == listed_tag:
                tail.reverse()
                return "".join(tail)
            else:
                tail.extend(current_tag)
                current_tag = []
            in_tag = False
        fhand.seek(-2, 1)

def xml_itemize(fhand, tag, num_items=1):
    '''It takes a xml file and it chunks it by the given key. It adds header if
    exists to each of the pieces. It is a generator'''
    fhand.seek(0, 2)
    end_file = fhand.tell()

    header = _get_xml_header(fhand, tag)
    tail   = _get_xml_tail(fhand, tag)
    section      = []
    current_tag  = []
    listed_tag_s = '<' + tag + '>'
    listed_tag_e = '</' + tag + '>'
    in_tag, in_section = False, False

    fhand.seek(0, 0)

    items_in_buffer = 0
    buffer = []
    while True:
        if end_file <= fhand.tell():
            break
        letter = fhand.read(1)
        if letter == '<':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        if in_section:
            section.append(letter)
        if letter == '>':
            if listed_tag_s == _remove_atributes_to_tag(current_tag):
                in_section = True
                section.extend(current_tag)
            elif listed_tag_e == _remove_atributes_to_tag(current_tag):
                items_in_buffer += 1
                buffer.extend(section)
                if items_in_buffer >= num_items:
                    yield  header + "".join(buffer) + tail
                    items_in_buffer = 0
                    buffer = []
                section = []
                in_section = False
            in_tag = False
            current_tag = []
    #is there any reamining buffer
    if buffer:
        yield  header + "".join(buffer) + tail
