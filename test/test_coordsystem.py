'''
Created on Oct 4, 2010

@author: peio
'''
from human_mut.genome import (CoordSystem, CodonPosition)
import unittest

class CodonPositionTest(unittest.TestCase):
    'It tests the CodonPosition class'

    def test_basic_codon(self):
        'It test the codon comparations'

        assert CodonPosition(0, 0) == CodonPosition(0, 0)
        assert CodonPosition(1, 1) == 4
        assert CodonPosition(0, 0) != CodonPosition(0, 1)
        assert CodonPosition(0, 0) < CodonPosition(0, 1)
        assert CodonPosition(0, 1) < CodonPosition(1, 0)
        assert CodonPosition(1, 0) > CodonPosition(0, 1)
        assert CodonPosition(1, 0) >= CodonPosition(0, 1)
        assert CodonPosition(1, 0) <= 3 <= CodonPosition(10, 1)

        assert CodonPosition(0, 1) + 2 == CodonPosition(1, 0)
        assert CodonPosition(0, 1) + 3 == CodonPosition(1, 1)
        assert 3 + CodonPosition(0, 1) == CodonPosition(1, 1)

        assert CodonPosition(-1, -1) == -4
        assert CodonPosition(0, 0) == 0
        assert 0 == CodonPosition(0, 0)
        assert not 0 != CodonPosition(0, 0)
        assert not CodonPosition(0, 0)
        assert CodonPosition(1, 1)- CodonPosition(0, 0) == 4

class CoordSystemTest2(unittest.TestCase):
    'It tests the coordinate transformations'

    def test_basic_coord(self):
        'It tests the basic coordinate transformations'

        #geno              111
        #        0123456789012

        #geno2   1234567890123

        #geno3   1111111
        #        6543210987654
        #
        #cdna      0123  4567
        #
        #prot        00  011
        #            01  201

        coord = CoordSystem(relations=[{'geno': [(2, 5), (8, 11)],
                                        'cdna': [(0, 3), (4, 7)]}])
        assert coord.transform(from_mol='cdna', to_mol='geno', position=6) == 10
        assert coord.transform(from_mol='cdna', to_mol='geno', position=0) == 2

        #more than two
        coord = CoordSystem(relations=[{'geno': [(0, 12)], 'geno2': [(1, 13)]},
                                       {'geno': [(2, 5), (8, 11)],
                                                     'cdna':[(0, 3), (4, 7)]},])

        assert coord.transform(from_mol='geno', to_mol='cdna', position=3) == 1
        assert coord.transform(from_mol='cdna', to_mol='geno', position=6) == 10
        assert coord.transform(from_mol='geno', to_mol='geno2',
                                                               position=1) == 2

        #with proteins
        coord = CoordSystem(relations=[{'geno': [(0, 12)], 'geno2': [(1, 13)]},
                                       {'geno': [(2, 5), (8, 11)],
                                                       'cdna':[(0, 3), (4, 7)]},
                                        {'cdna':[(2, 6)],
                         'prot':[(CodonPosition(0, 0), CodonPosition(1, 1))]}])

        assert coord.transform(from_mol='geno', to_mol='cdna', position=3) == 1
        assert coord.transform(from_mol='cdna', to_mol='geno', position=6) == 10
        assert coord.transform(from_mol='geno', to_mol='geno2',
                                                              position=11) == 12
        assert coord.transform(from_mol='geno', to_mol='prot',
                                                              position=5) == 1
        assert coord.transform(from_mol='prot', to_mol='geno',
                                                              position=2) == 8

        #geno              111
        #        0123456789012

        #geno3   1111111
        #        6543210987654

        #reversed
        coord = CoordSystem(relations=[{'geno': [(0, 12)], 'geno3': [(16, 4)]}])
        assert coord.transform(from_mol='geno', to_mol='geno3',
                                                             position=2) == 14

        #geno              111
        #        0123456789012

        #cdna     876  543210
        #prot      54  3210

        coord = CoordSystem(relations=[{'geno': [(1, 3), (6, 11)],
                                                       'cdna':[(8, 6), (5, 0)]},
                                        {'cdna':[(2, 7)],
                         'prot':[(CodonPosition(0, 0), CodonPosition(1, 2))]}])
        assert coord.transform(from_mol='geno', to_mol='prot',
                                                             position=2) == 5


    def test_transformation(self):
        'It test the transformation between fragments'
        coord = CoordSystem(relations=[])

        assert coord._transform_pos((1, 10), (2, 11), 2) == 3
        assert coord._transform_pos((1, 10), (11, 2), 2) == 10
        assert coord._transform_pos((0, 12), (16, 4), 2) == 14

        #if we deal with CodonPosition we get the right answer and type
        result = coord._transform_pos((2, 7), (CodonPosition(2, 1),
                                               CodonPosition(4, 0)), 4)
        expected = CodonPosition(3, 0)
        assert result == expected
        assert type(result) == type(expected)

        result = coord._transform_pos((CodonPosition(2, 1),
                                       CodonPosition(4, 0)),
                                       (2, 7), CodonPosition(3, 1))
        expected = 5
        assert result == expected
        assert type(result) == type(expected)

        result = coord._transform_pos((CodonPosition(4, 0),
                                       CodonPosition(2, 1)),
                                       (2, 7), CodonPosition(3, 2))
        expected = 3
        assert result == expected
        assert type(result) == type(expected)


    def test_merge_overlaping_segments(self):
        'It test the segment merging'
        coord = CoordSystem(relations=[])
        limits = [(0, 2)]
        expected = [(0, 2)]
        assert coord._merge_overlaping_segments(limits) == expected

        limits = [(0, 2), (4, 6)]
        expected = [(0, 2), (4, 6)]
        assert coord._merge_overlaping_segments(limits) == expected

        limits = [(0, 2), (4, 6), (0, 8)]
        expected = [(0, 8)]
        assert coord._merge_overlaping_segments(limits) == expected

    def test_intersec_segments(self):
        'It test the segment intersection'
        coord = CoordSystem(relations=[])

        limits1 = [(0, 2)]
        limits2 = [(2, 8)]
        expected = [(0, 1), (2, 2), (3, 8)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(1, 10)]
        limits2 = [(5, 6)]
        expected = [(1, 4), (5, 6), (7, 10)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(1, 10)]
        limits2 = [(1, 10)]
        expected = [(1, 10)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(1, 10)]
        limits2 = [(2, 9)]
        expected = [(1, 1), (2, 9), (10, 10)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(4, 12)]
        limits2 = [(6, 8)]
        expected = [(4, 5), (6, 8), (9, 12)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(4, 8)]
        limits2 = [(6, 12)]
        expected = [(4, 5), (6, 8), (9, 12)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(0, 8)]
        limits2 = [(0, 12)]
        expected = [(0, 8), (9, 12)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(0, 8)]
        limits2 = [(2, 8)]
        expected = [(0, 1), (2, 8)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(0, 2)]
        limits2 = [(3, 8)]
        expected = [(0, 2), (3, 8)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(0, 2)]
        limits2 = [(4, 8)]
        expected = [(0, 2), (4, 8)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(0, 12)]
        limits2 = [(2, 5), (8, 11)]
        expected = [(0, 1), (2, 5), (6, 7), (8, 11), (12, 12)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(12, 0)]
        limits2 = [(12, 0)]
        expected = [(12, 0)]
        assert coord._intersec_segments(limits1, limits2) == expected

        limits1 = [(12, 6), (5, 0)]
        limits2 = [(13, 7)]
        expected = [(13, 13), (12, 7), (6, 6), (5, 0)]
        assert coord._intersec_segments(limits1, limits2) == expected

        #cdna     876  543210
        #prot      54  3210
        limits1 = [(8, 6), (5, 0)]
        limits2 = [(7, 2)]
        expected = [(8, 8), (7, 6), (5, 2), (1, 0)]
        assert coord._intersec_segments(limits1, limits2) == expected


if __name__ == "__main__":
#    test = 'CoordSystemTest2.test_intersec_segments'
#    test = 'CoordSystemTest2.test_merge_overlaping_segments'
    #import sys;sys.argv = ['', test]

    unittest.main()
