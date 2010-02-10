'''Tests for the collections module'''

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

import unittest, StringIO
from biolib.utils.collections_ import (item_context_iter,
                                       RequiredPosition, list_pairs_iter)

class _Locatable():
    'A class for locatable items'
    def __init__(self, location, reference):
        'The init only requires a location'
        self.location = location
        self.reference = reference
    def __str__(self):
        'Show myself'
        return '%s %d' % (self.reference, float(self.location))
    def __repr__(self):
        'same as str'
        return self.__str__()

def _check_item_context_iter(item_context_iter_, expected):
    'It check that the given iter match the expected items'
    for index, (item, context) in enumerate(item_context_iter_):
        assert expected[index][0][0] == item.reference
        assert expected[index][0][1] == item.location
        assert len(context) == len(expected[index][1])
        for context_index, (reference, location) in \
                                            enumerate(expected[index][1]):
            assert context[context_index].reference == reference
            assert context[context_index].location == location

class ItemContextIterTest(unittest.TestCase):
    'It checks the filtering methods.'

    @staticmethod
    def test_item_context_iter():
        'It tests the item contex iterator with a window'
        #a starndard one
        items = [_Locatable(1, 'a'), _Locatable(2, 'a'), _Locatable(4, 'a'),
                 _Locatable(6, 'a'), _Locatable(1, 'b'), _Locatable(5, 'b')]
        items = iter(items)
        expected = [(('a', 1), (('a', 1), ('a', 2))),
                    (('a', 2), (('a', 1), ('a', 2))),
                    (('a', 4), (('a', 4),)),
                    (('a', 6), (('a', 6),)),
                    (('b', 1), (('b', 1),)),
                    (('b', 5), (('b', 5),)),]
        _check_item_context_iter(item_context_iter(items, window=3), expected)

        #an empty one
        items = []
        items = iter(items)
        expected = []
        _check_item_context_iter(item_context_iter(items, window=3), expected)

    @staticmethod
    def test_no_window():
        'If no window is given all items from a reference will be returned'
        items = [_Locatable(1, 'a'), _Locatable(5, 'b')]
        items = iter(items)
        expected = [(('a', 1), (('a', 1),)),
                    (('b', 5), (('b', 5),)),]
        _check_item_context_iter(item_context_iter(items), expected)

        items = [_Locatable(1, 'a'), _Locatable(2, 'a'), _Locatable(6, 'a'),
                 _Locatable(1, 'b'), _Locatable(5, 'b')]
        items = iter(items)
        expected = [(('a', 1), (('a', 1), ('a', 2), ('a', 6))),
                    (('a', 2), (('a', 1), ('a', 2), ('a', 6))),
                    (('a', 6), (('a', 1), ('a', 2), ('a', 6))),
                    (('b', 1), (('b', 1), ('b', 5))),
                    (('b', 5), (('b', 1), ('b', 5))),]
        _check_item_context_iter(item_context_iter(items), expected)

REQUIRED_POS = '''
aaa    2
aaa    3
aaa    5
aaa    8
bbb    2
bbb    7'''

REQUIRED_POS2 = '''
aaa    3
aaa    2
aaa    8
bbb    2
bbb    7'''

class RequiredPositionsTest(unittest.TestCase):
    'It checks the FileCachedRequiredPositions methods'
    @staticmethod
    def test_required_pos():
        'Basic class usage'
        fhand = StringIO.StringIO()
        fhand.write(REQUIRED_POS)
        fhand.flush()
        req_pos = RequiredPosition(fhand)
        assert req_pos['aaa', '2']
        assert req_pos['aaa', '3']
        assert not req_pos['aaa', '4']
        assert not req_pos['aaa', '4']
        assert req_pos['aaa', '5']
        assert not req_pos['aaa', '4']
        assert req_pos['bbb', '2']

    @staticmethod
    def test_required_pos_no_sorted():
        ''' It test if the file is not ordered  or if you ask in an unsorted
        way'''
        # File unordered
        fhand = StringIO.StringIO()
        fhand.write(REQUIRED_POS2)
        fhand.flush()
        req_pos = RequiredPosition(fhand)
        assert not req_pos['aaa', '2']

        # you ask unordered
        fhand = StringIO.StringIO()
        fhand.write(REQUIRED_POS)
        fhand.flush()
        req_pos = RequiredPosition(fhand)
        assert req_pos['aaa', '2']
        assert req_pos['aaa', '2']

    @staticmethod
    def test_pairs():
        'We can get all pair for the elements of a list'
        items = [1, 2, 3, 4]
        pairs = list(list_pairs_iter(items))
        assert list(pairs) == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

if __name__ == "__main__":
    unittest.main()