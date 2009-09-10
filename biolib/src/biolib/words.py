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
from biolib.collections_ import SparseVector
from biolib.rbtree import RBDict

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
