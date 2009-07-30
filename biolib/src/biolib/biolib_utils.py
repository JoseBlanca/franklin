'''
Created on 2009 api 30

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


import tempfile, shutil
import os, re, math

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



def split_long_sequences(seq_iter, max_length):
    'It splits thequences in this iterator taht excedes the max_length permited'
    for seq in seq_iter:
        seq_len = len(seq)
        if seq_len > max_length:
            seqs = _split_seq(seq, max_length)
            for seq in seqs:
                yield seq
        else:
            yield seq

def _split_seq(seq, maxlength):
    'It split sequences in a '
    seq_len       = len(seq)
    num_sequences = seq_len / maxlength

    if  seq_len % maxlength != 0 and num_sequences == 1:
        num_sequences += 1

    # Calculate start and end of the new sequences with the original sequence
    # coordinates
    start_ends = _calculate_divisions(seq_len, num_sequences)
    for i, (start, end) in enumerate(start_ends):
        if seq.qual is not None:
            qual = seq.qual[start:end]
        else:
            qual = None
        yield seq.copy(seq=seq.seq[start:end], qual=qual,
                        name='%s_%d' % (seq.name, i+1))

def _calculate_divisions(length, splits):
    '''It calculates the length of each new seq.
    It returns the start and end of the new sequences. giving the coordinates
    of the original sequence.
    '''
    num_fragments1 = length % splits
    num_fragments2 = splits - num_fragments1
    new_length2 = length // splits
    new_length1 = new_length2 + 1
    res = ((num_fragments1, new_length1), (num_fragments2, new_length2))

    lengths = []
    for num_fragments, length in res:
        for i in range(num_fragments):
            lengths.append(length)
    # Now I calculate the start and end of each sequence
    r_length   = 0
    start_ends = []
    for i, length in enumerate(lengths):
        if i == 0:
            start_ends.append((0, length))
        else:
            start_ends.append((r_length, r_length + length))
        r_length += length
    return start_ends




def get_start_end(location):
    '''It accepts an int, Location or tuple and it returns the start, end,
    forward and strand.'''
    #int
    if isinstance(location, int):
        start = location
        end = location
    #tuple
    elif isinstance(location, tuple):
        start = location[0]
        end = location[1]
    #location
    else:
        start = location.start
        end = location.end
    return start, end



def remove_from_orf(orf_dna, orf_prot, aminoacid='X'):
    ''' It removes an aminoaacid from dna and protein seq'''
    dna  = []
    prot = []
    pos       = 0
    for letter in orf_prot:
        if letter.upper() != aminoacid:
            prot.append(letter)
            dna.append(orf_dna[pos:pos + 3])
        pos += 3
    return "".join(dna), "".join(prot)

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

def xml_itemize(fhand, tag):
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
                yield  header + "".join(section) + tail
                section = []
                in_section = False
            in_tag = False
            current_tag = []

def get_safe_fname(directory, prefix, suffix):
    '''It looks if the name exits in this directory and adds a number to the end
    to be a threath same name'''

    if not os.path.exists('%s/%s.%s' % (directory, prefix, suffix)):
        return '%s/%s.%s' % (directory, prefix, suffix)
    for i in range(10000):
        if not os.path.exists('%s/%s.%d.%s' %(directory, prefix, i,
                                                     suffix)):
            return '%s/%s.%d.%s' % (directory, prefix, i, suffix)

    msg  = 'There are more than 100000 outputs of this kind is this directory'
    raise ValueError(msg)

class _FileItemGetter(object):
    '''Syntax sugar for the File index.

    It clarifies the use of the FileIndex for every item type.
    '''
    def __init__(self, fhand, items):
        '''It requires an fhand and a dict of items.

        In each item it should be the item start in the file and the
        length in bytes
        '''
        self._fhand = fhand
        self._items = items

    def __getitem__(self, key):
        'It returns one item from the file'
        start, length = self._items[key]
        self._fhand.seek(start)
        return self._fhand.read(length)

    def items(self):
        'it yields the item keys'
        for item in self._items:
            yield item

class FileIndex(object):
    '''This class is used to index the items present in a file.

    We can consider a file as a stream of items separated by certain patterns.
    This class assigns a key and a kind to each of these items and it creates
    an interface to get these items in a convenient way. It works like:
    index = FileIndex(fhand, item_start_patterns=['pattern1', 'pattern2'])
    index['default_type'][key]
    '''
    _default_item = 'item'
    def __init__(self, fhand, item_start_patterns, key_patterns,
                 type_patterns=None):
        '''It indexes the given fhand.

        The items in the file are separated by one of the patterns given in the
        item_start_patterns list.
        The key for each pattern is found with the regex present in the
        key_patterns list. Optionally a dict with type and patterns for each
        type can be given. If these types are not given all the items will have
        the type 'item'.
        The subitems patterns can't be regex, they're just chars
        '''
        self._fhand = fhand
        self._item_start_patterns = item_start_patterns
        self._key_patterns = key_patterns
        self._type_patterns = type_patterns
        self._subitem_patterns = [os.linesep]
        self._all_patterns_to_regex()
        self._index = None
        self._create_index()
        self._getters = None
        self._create_properties()

    @staticmethod
    def _patterns_to_regex(patterns):
        'Given a list of patterns (str or regex) it returns a list of regex'
        new_patterns = []
        for pattern in patterns:
            if 'match' not in dir(pattern): #if it isn't a regex
                pattern = re.compile(pattern)
            new_patterns.append(pattern)
        return new_patterns

    def _all_patterns_to_regex(self):
        '''It transforms all patterns in regular expressions.

        They could be either str or regex.
        '''
        patterns_to_regex = self._patterns_to_regex
        self._item_start_patterns = patterns_to_regex(self._item_start_patterns)
        self._key_patterns = patterns_to_regex(self._key_patterns)
        #the subitem patterns a not transformed because they're just chars

        #the type patterns can be a list, dict or None
        type_patterns = self._type_patterns
        if isinstance(type_patterns, dict):
            if type_patterns is not None:
                for type_ in type_patterns:
                    type_patterns[type_] = \
                                         patterns_to_regex(type_patterns[type_])
        elif isinstance(type_patterns, list):
            type_patterns = patterns_to_regex(type_patterns)
        self._type_patterns = type_patterns

    @staticmethod
    def _patterns_match_in_string(string, patterns):
        '''It returns True if a match is found in the string for any of the
        given patterns'''
        for pattern in patterns:
            if pattern.match(string):
                return True
        return False

    @staticmethod
    def _subitems_in_file(fhand, patterns):
        '''It yields all subitems in a file.

        The subitems will be divided by the given patterns.
        '''
        current_subitem = ''
        while True:
            char = fhand.read(1)
            if not char:
                #file done
                yield current_subitem
                break
            current_subitem += char
            if char in patterns:
                yield current_subitem
                current_subitem = ''

    def _items_in_fhand(self):
        '''It yields all the items in the fhand.

        It yields tuples with (key, type, start, end)
        '''
        #a file is composed by items, each item is composed by sub items.
        #for instance a fasta file is build by sequence items, each of them is
        #build by lines. An xml file is build by items and each of them is
        #composed by <subitems>
        fhand = self._fhand
        subitem_patterns = self._subitem_patterns
        item_start_patterns = self._item_start_patterns
        key_patterns = self._key_patterns
        type_patterns = self._type_patterns
        patterns_match_in_string = self._patterns_match_in_string
        subitems_in_file = self._subitems_in_file

        def update_key(key, string):
            '''If the string match in any of the key_patterns it returns the key

            If the key was not None it raises an error because two keys have
            been found in the same item.
            '''
            for pattern in key_patterns:
                match = pattern.search(string)
                if match:
                    new_key = match.group(1)
                    if key is None:
                        return new_key
                    else:
                        raise RuntimeError('Two keys found in an item: %s,%s' %
                                           (key, new_key))
            return key

        def update_type(type_, string):
            '''If the string match in any of the type_patterns it returns the
            type

            If the type_ was not None it raises an error because two types have
            been found in the same item.
            '''
            #type_patterns can be None, dict or list
            new_type = None
            if type_patterns is None:
                return self._default_item
            elif isinstance(type_patterns, dict):
                for known_type, patterns in type_patterns.items():
                    for pattern in patterns:
                        match = pattern.search(string)
                        if match:
                            new_type = known_type
                            break
            elif isinstance(type_patterns, list):
                for pattern in type_patterns:
                    match = pattern.search(string)
                    if match:
                        new_type = match.group(1)
            if new_type is None and type_ is not None:
                return type_
            elif new_type is not None and type_ is None:
                return new_type
            elif new_type is not None and type_ is not None:
                msg = 'Two types found in an item: %s,%s' % (type_, new_type)
                raise RuntimeError(msg)
            return type_

        fhand.seek(0)
        current_item_start = 0
        previous_position = 0
        current_key = None
        current_type = None
        for subitem in subitems_in_file(fhand, subitem_patterns):
            if previous_position and patterns_match_in_string(subitem,
                                                           item_start_patterns):
                length = previous_position - current_item_start
                yield current_item_start, length, current_type, current_key
                current_item_start = previous_position
                current_key = None
                current_type = None
            current_key = update_key(current_key, subitem)
            current_type = update_type(current_type, subitem)
            previous_position = fhand.tell()
        else:
            length = previous_position - current_item_start
            yield current_item_start, length, current_type, current_key

    def _create_index(self):
        'It creates the index for the file'
        index = {}
        for item in self._items_in_fhand():
            start, length, type_, key = item
            if type_ not in index:
                index[type_] = {}
            index[type_][key] = (start, length)
        self._index = index

    def _create_properties(self):
        'It creates a property for every type found in the file'
        getters = {}
        for type_, items in self._index.items():
            getters[type_] = _FileItemGetter(self._fhand, items)
        self._getters = getters

    def __getitem__(self, type_):
        '''It returns file items'''
        if self._type_patterns is None:
            return self._getters[self._default_item][type_]
        else:
            return self._getters[type_]

    def types(self):
        'It yields all the type names'
        for type_ in self._index:
            yield type_

class FileCachedList(object):
    '''A list cached in a file.

    The main aim is to store really big lists in files with an iterator
    interface.
    '''
    def __init__(self, type_):
        '''It creates a FileCachedList.

        It requires the type of objects that will be stored.
        '''
        self._type = type_
        #the file to store the data
        self._wfhand = tempfile.NamedTemporaryFile()
        #the file to read the data
        self._rfhand = None

    def append(self, item):
        'It writes one element at the file end'
        self._wfhand.write('%s\n' % str(item))

    def extend(self, items):
        'It adds a bunch of items from a list or from a FileCachedList'
        if 'items' in dir(items):
            items = items.items()
        for item in items:
            self.append(item)

    def items(self):
        'It yields all items'
        self._wfhand.flush()
        rfhand = open(self._wfhand.name)
        for line in rfhand:
            if len(line) == 0:
                raise StopIteration
            yield self._type(line.strip())
