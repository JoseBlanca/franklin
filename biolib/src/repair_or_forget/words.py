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

'''Utilities to deal with sequence words

Created on 03/09/2009

@author: jose
'''
import array
from biolib.rbtree import RBDict
import pysparse

def cluster_seqs_by_words(seqs, wsize, count=False):
    '''It looks for sequences that share at least one word'''
    code = {'A':'0', 'T':'1', 'C':'2', 'G':'3'}
    # We need
    nwords = int('3'* wsize, 4) + 1
    word_index = SparseVector(nwords, store_non_int=not count)
    seq_names = []

    for seq_number, seq in enumerate(seqs):
        print seq_number
        seq_names.append(seq.name)
        seq   = [code.get(letter, 'x') for letter in str(seq.seq)]
        if count:
            analyze_words_count(seq, wsize, word_index)
        else:
            analyze_words(seq, wsize, word_index, seq_number)
    return word_index, seq_names

def cluster_seqs_by_words2(seqs, wsize, count=False):
    '''It looks for sequences that share at least one word'''
    code = {'A':'0', 'T':'1', 'C':'2', 'G':'3'}
    # We need
    word_index = RBDict()

    seq_names = []

    for seq_number, seq in enumerate(seqs):
        print seq_number
        seq_names.append(seq.name)
        seq   = [code.get(letter, 'x') for letter in str(seq.seq)]
        if count:
            analyze_words_count(seq, wsize, word_index)
        else:
            analyze_words(seq, wsize, word_index, seq_number)
    return word_index, seq_names

def analyze_words_count(seq, wsize, word_index):
    '''It count repeats of each word and saves the
     information  in word_index, so it does not need to return it.

     We store the words codified in base 10 to be able to use the list index
     to get/set them
     '''
    for word in words_in_seq(seq, wsize):
        word = int("".join(word), 4)
        word_index[word] += 1

def analyze_words(seq, wsize, word_index, seq_number):
    '''it group sequences if they share a word.

    It stores in word_index the word codified as the index.
    The content of each index is an array object holding the list of sequences
    that share this word'''
    for word in words_in_seq(seq, wsize):
        word = int("".join(word), 4)
        if word_index[word] is None:
            word_index[word] = array.array('I')

        this_array = word_index[word]
        if seq_number not in this_array:
            this_array.append(seq_number)

def words_in_seq(seq, wsize):
    '''It yields all the words in a sequuence of a given size, that are not masked.
     Before using this function you must codificate the masked nucleotides
     as 'x'  '''
    for i in range(len(seq) - wsize + 1):
        word = seq[i:i+wsize]
        if 'x' not in word:
            yield word

BASE2 = "01"
BASE4 = '0123'
BASE10 = "0123456789"
BASE16 = "0123456789ABCDEF"
BASE62 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"

def baseconvert(number,fromdigits,todigits):
    """ converts a "number" between two bases of arbitrary digits

    The input number is assumed to be a string of digits from the
    fromdigits string (which is in order of smallest to largest
    digit). The return value is a string of elements from todigits
    (ordered in the same way). The input and output bases are
    determined from the lengths of the digit strings. Negative
    signs are passed through.

    decimal to binary
    >>> baseconvert(555,BASE10,BASE2)
    '1000101011'

    binary to decimal
    >>> baseconvert('1000101011',BASE2,BASE10)
    '555'

    integer interpreted as binary and converted to decimal (!)
    >>> baseconvert(1000101011,BASE2,BASE10)
    '555'

    base10 to base4
    >>> baseconvert(99,BASE10,"0123")
    '1203'

    base4 to base5 (with alphabetic digits)
    >>> baseconvert(1203,"0123","abcde")
    'dee'

    base5, alpha digits back to base 10
    >>> baseconvert('dee',"abcde",BASE10)
    '99'

    decimal to a base that uses A-Z0-9a-z for its digits
    >>> baseconvert(257938572394L,BASE10,BASE62)
    'E78Lxik'

    ..convert back
    >>> baseconvert('E78Lxik',BASE62,BASE10)
    '257938572394'

    binary to a base with words for digits (the function cannot convert this back)
    >>> baseconvert('1101',BASE2,('Zero','One'))
    'OneOneZeroOne'

    """

    if str(number)[0]=='-':
        number = str(number)[1:]
        neg=1
    else:
        neg=0

    # make an integer out of the number
    x=long(0)
    for digit in str(number):
       x = x*len(fromdigits) + fromdigits.index(digit)

    # create the result in base 'len(todigits)'
    res=""
    while x>0:
        digit = x % len(todigits)
        res = todigits[digit] + res
        x /= len(todigits)
    if neg:
        res = "-"+res
    return res

def filter_low_abundant_words(word_abundance, min_abundance, wsize):
    filtered_words = []
    letters = ['A', 'T', 'C', 'G']
    for word, abundance in word_abundance.non_empty_items():
        if abundance > min_abundance:
            word = str(baseconvert(word, BASE10, BASE4)).rjust(wsize, '0')
            word = [letters[int(numb)] for numb in word]
            word = ''.join(word)
            filtered_words.append((word, abundance))
    def cmp(abundance1, abundance2):
        return int(abundance1[1] - abundance2[1])
    filtered_words = sorted(filtered_words, cmp)
    return filtered_words


class SparseVector(object):
    'A sparse vector collection to hold nearly empty vectors.'
    def __init__(self, nelements, size_hint=None, store_non_int=False):
        '''It initializes the vector.

        We require the number of elements that the vector should hold.
        size_hint is the number of estimated elements that the vector should
        hold at the end. If given expensive memory reallocation will be saved.
        '''
        self._size = nelements
        if size_hint:
            self._spv = pysparse.spmatrix.ll_mat(nelements, 1, size_hint)
        else:
            self._spv = pysparse.spmatrix.ll_mat(nelements, 1)
        if store_non_int:
            self._cargo = []
            self._cargo.append(None)   #we dont' use the 0 as cargo index
        else:
            self._cargo = None
        self._position = 0

    def __setitem__(self, index, value):
        'It sets one item in the vector'
        cargo = self._cargo
        if cargo is None:
            self._spv[index, 0] = value
        else:
            #was an item already stored in that position?
            cargo_index = self._spv[index, 0]
            if cargo_index:
                self._cargo[cargo_index] = value
            else:
                self._cargo.append(value)
                cargo_index = len(self._cargo) - 1
                self._spv[index, 0] = cargo_index

    def __getitem__(self, index):
        'It returns one element'
        cargo = self._cargo
        if cargo is None:
            return self._spv[index, 0]
        else:
            return cargo[int(self._spv[index, 0])]
    def __iter__(self):
        'The iterable protocol'
        self._position = 0
        return self

    def next(self):
        'It returns the next item'
        position = self._position
        if position < self._size:
            value = self[self._position]
            self._position += 1
            return value
        else:
            raise StopIteration

    def __len__(self):
        'It returns the number of elements'
        return self._size

    def __str__(self):
        'It prints the vector'
        toprint = '['
        for item in self:
            toprint += str(item)
            toprint += ', '
        toprint += ']'
        return toprint

    def non_empty_items(self):
        'It return a iterator with non empty items'
        keys = self._spv.keys()[0]
        last_item = len(keys) - 1
        for index, key in enumerate(keys):
            if self._cargo:
                yield key, self._cargo[int(self._spv[key, 0])]
            else:
                yield key, self._spv[key, 0]
            if index == last_item:
                raise StopIteration

    def non_empty_values(self):
        'It return a iterator with non empty items'
        for item in self.non_empty_items():
            yield item[1]
