'''It checks that the contig code is ok.'''
# disable: Too many lines in module
# pylint: disable-msg=C0302

import unittest
from biolib.contig import Location, NonStaticParentLocation, \
                          LocatableSequence, _reparent_location, \
                          locate_sequence, Contig
from test.test_utils import Seq, SeqWithQuality, SeqRecord, Seqmut

# pylint: disable-msg=R0201
# disable to many public methods
# pylint: disable-msg=R0904
# disable to many statements
# pylint: disable-msg=R0915
# disable: Except doesn't do anything
# pylint: disable-msg=W0704

class LocationTests(unittest.TestCase):
    '''It tests all the functionallity offered by the Location class.'''
    def test_basic_functionality(self):
        '''It tests the start and end properties and the length'''
        #with a parent sequence start and end are optional.
        seq = '0123456789'
        loc = Location(parent=seq)
        assert loc == Location(start=0, end=len(seq) - 1, parent=seq)

        #we can give start and/or end
        loc = Location(start=1, parent=seq)
        assert loc == Location(start=1, end=len(seq) - 1, parent=seq)

        loc = Location(end=5, parent=seq)
        assert loc == Location(start=0, end=5, parent=seq)

        #the parent is optional
        loc = Location(end=5)
        assert loc == Location(start=0, end=5)

        #the start should be lower than the end
        self.failUnlessRaises(ValueError, Location, start=1, end=0)
        self.failUnlessRaises(ValueError, Location, start= - 1, end=10)
        self.failUnlessRaises(ValueError, Location, start=10, end=20,
                              strand=2)

        #some erros should be raised
        self.failUnlessRaises(ValueError, Location, start=1, end=0)
        self.failUnlessRaises(ValueError, Location, start= - 1, end=10)
        self.failUnlessRaises(ValueError, Location, start=10, end=20,
                                                                  strand=2)

    def test_len(self):
        '''It checks that the length method works as expected'''
        loc = Location(end=5)
        assert len(loc) == 6

        seq = '0123456789'
        #it behaves as a sequence
        loc = Location(parent=seq)
        assert len(loc) == len(seq)

        #when the location does not starts at 0
        loc = Location(start=2, end=5)
        assert len(loc) == 4

    def test_reverse(self):
        '''Test that the Location is not mutable, it returns a new instance.
        '''
        # 012345
        #  |---|
        loc = Location(start=1, end=5, strand=1)
        loc2 = loc[::-1]
        assert loc2 == Location(start=1, end=5, strand=1, forward=False)

        #reverse with a parent sequence of 10 residues
        seq = '0123456789'
        # 0123456789
        #  |---|
        #     |---|
        loc = Location(start=1, end=5, strand=-1, forward=False, parent=seq)
        loc2 = loc[::-1]
        assert loc2 == Location(start=1, end=5, strand=-1, parent=seq)

        # 01234567890
        #      |---|
        #  |---|
        seq = '01234567890'
        loc = Location(start=5, end=9, parent=seq)
        loc2 = loc[::-1]
        assert loc2 == Location(start=5, end=9, forward=False, parent=seq)

        try:
            #the lower than -1 strides are not supported
            loc2 = loc[::-2]
            self.fail('IndexError expected')
        except IndexError:
            pass

    def test_complement(self):
        '''Test that the Location is not mutable, when complemented it returns
        a new instance.'''
        seq = '012345'
        loc = Location(start=1, end=5, strand=1, parent=seq)
        loc2 = loc.complement()
        assert loc2 == Location(start=1, end=5, strand=-1, parent=seq)

        loc = Location(start=1, end=5, strand= - 1, forward=False)
        loc2 = loc.complement()
        assert loc2 == Location(start=1, end=5, strand=1, forward=False)

        loc = Location(start=1, end=5, parent=seq)
        loc2 = loc.complement()
        assert loc2 == Location(start=1, end=5, parent=seq)

    def test_slices(self):
        '''We can get a sliced location, in a similar way as if it was a
        list. The parent is not affected by the slicing.'''
        #a standard slice inside the location
        #54321 <- negative index
        #01234
        #  ---
        #  +++*
        #  012
        loc = Location(start=2, end=4, strand=1)
        loc2 = loc[2:5]
        assert loc2 == Location(start=2, end=4, strand=1)
        #with negative indexes
        loc = Location(start=2, end=4, strand=1, forward=False)
        loc2 = loc[-3:]
        assert loc2 == Location(start=2, end=4, strand=1, forward=False)

        #only partially inside the location
        #54321 <- negative index
        #01234567
        #  ---
        # +++*
        # 012
        loc = Location(start=2, end=4, strand=1)
        loc2 = loc[1:4]
        assert loc2 == Location(start=2, end=3, strand=1)
        #with negative indexes
        loc = Location(start=2, end=4, strand=1)
        loc2 = loc[-4:-1]
        assert loc2 == Location(start=2, end=3, strand=1)

        #01234567
        #  ---
        #   +++*
        #   01
        loc = Location(start=2, end=4, strand=1)
        loc2 = loc[3:6]
        assert loc2 == Location(start=3, end=4, strand=1)

        #01234567
        #  ---
        # +++++*
        # 0123
        loc = Location(start=2, end=4, strand=1)
        loc2 = loc[1:6]
        assert loc2 == Location(start=2, end=4, strand=1)

        #an edge case
        #01234567
        #  ---
        #    +++*
        #    0
        loc = Location(start=2, end=4, strand=1)
        loc2 = loc[4:7]
        assert loc2 == Location(start=4, end=4, strand=1)

        #an edge case
        #54321 <- negative index
        #01234
        #  ---
        #++++*
        #0123
        loc = Location(start=2, end=4, strand=1)
        loc2 = loc[:4]
        assert loc2 == Location(start=2, end=3, strand=1)
        #with negative index
        loc2 = loc[-5:-1]
        assert loc2 == Location(start=2, end=3, strand=1)
        
        #a single int as index
        #54321 <- negative index
        #01234
        #  ---
        #   ^
        #   0
        loc = Location(start=2, end=4, strand=1)
        loc2 = loc[3]
        assert loc2 == Location(start=3, end=3, strand=1)
        #with negative index
        loc2 = loc[-2]
        assert loc2 == Location(start=3, end=3, strand=1)

        #None returns everything
        loc = Location(start=10, end=20, strand=1)
        loc2 = loc[:None]
        assert loc2 == loc
        loc = Location(start=10, end=20, strand=1)
        loc2 = loc[:]
        assert loc2 == loc

        #a stride different than 1
        #01234567890
        #-----------
        # +   + *
        # 01234
        loc = Location(start=0, end=10, strand=1)
        try:
            loc2 = loc[1:7:4]
            self.fail('IndexError expected')
        except IndexError:
            pass

        #negative index with a lengthier parent sequence
        #a parent sequence of 10 residues
        seq = '0123456789'
        # 0987654321 <- negative index
        # 0123456789
        #  -----
        #    -----*
        #    01234 <- new reference
        loc = Location(start=1, end=5, strand=-1, parent=seq)
        loc2 = loc[3:8]
        assert loc2 == Location(start=3, end=5, strand=-1, parent=seq)
        #with the negative index
        loc2 = loc[-7:-2]
        assert loc2 == Location(start=3, end=5, strand=-1, parent=seq)

        #index errors
        #outside the location
        loc = Location(start=10, end=20, strand=1)
        try:
            loc2 = loc[5:10]
            self.fail('IndexError expected')
        except IndexError:
            pass
        try:
            loc2 = loc[21:27]
            self.fail('IndexError expected')
        except IndexError:
            pass
        try:
            #a list would return an empty list,
            #but a biolocation is not mutable and it has not sense
            #to have an empty location
            loc2 = loc[5:5]
            self.fail('IndexError expected')
        except IndexError:
            pass
        try:
            #outside with an int
            loc2 = loc[27]
            self.fail('IndexError expected')
        except IndexError:
            pass
        try:
            #the -1 stride is not supported with start and stop
            loc2 = loc[10:20:-1]
            self.fail('IndexError expected')
        except IndexError:
            pass

    def test_apply_to_parent(self):
        '''The Location can be used to transform the parent sequence.'''
        #It should work with static and non-static parents
        seq = '0123456789'
        #static parents
        loc = Location(start=2, end=6, forward=True, parent=seq)
        assert loc.parent is seq
        assert loc.apply_to_parent() == seq[2:7]
        #with reverse
        loc = Location(start=2, end=6, forward=False, parent=seq)
        assert loc.apply_to_parent() == seq[2:7][::-1]

        #the static parent keeps the parent and can apply the transformation
        #after being sliced
        loc = Location(parent=seq)
        loc2 = loc[2:7][::-1]
        assert loc2.apply_to_parent() == seq[2:7][::-1]

    def test_coordinate_transformations(self):
        '''It transforms from the Location to the parent coordinates.'''
        # 0123456789
        #   01234
        seq = '0123456789'
        loc = Location(start=2, end=6, parent=seq)
        coord_pairs = ((1, 3), (slice(1, 4), slice(3, 6)),
                       (slice(None, None), slice(None, None)))
        for coord_pair in coord_pairs:
            loc_coord = coord_pair[0]
            parent_coord = coord_pair[1]
            assert loc.get_parent_index(loc_coord) == parent_coord
            assert loc.get_location_index(parent_coord) == loc_coord
        
        #some indexes outside the Location
        # 0123456789
        #   0123
        loc = Location(start=2, end=5, parent=seq)
        assert loc.get_location_index(slice(1, 5)) == slice(None, 3)
        assert loc.get_location_index(slice(2, 6)) == slice(0, None)
        assert loc.get_location_index(1) is None
        assert loc.get_location_index(6) is None
        #now oustide the location but with a reversed one
        # 0123456789
        #   3210
        loc = Location(start=2, end=5, forward=False, parent=seq)
        assert loc.get_location_index(slice(1, 5)) == slice(0, None)
        assert loc.get_location_index(slice(2, 6)) == slice(None, 3)
        assert loc.get_location_index(1) is None
        assert loc.get_location_index(6) is None

    def test_overlaps(self):
        '''Locations overlaps'''
        loc_pairs = (((2, 5), (1, 3), True), ((4, 5), (1, 6), True),
                     ((1, 5), (2, 6), True), ((1, 5), (2, 4), True),
                     ((0, 5), (6, 10), False), ((10, 55), (6, 8), False))
        for loc_pair in loc_pairs:
            loc1 = Location(start=loc_pair[0][0], end=loc_pair[0][1])
            loc2 = Location(start=loc_pair[1][0], end=loc_pair[1][1])
            tuple2 = loc_pair[1]
            #it works with locations
            assert loc1.overlaps(loc2) == loc_pair[2]
            #and with tuples (start, end)
            assert loc1.overlaps(tuple2) == loc_pair[2]

        #this test was added to catch an obscure error in a snp processing
        #tool
        loc1 = Location(start=1, end=303)
        assert not loc1.overlaps((305, None))

    def test_contains(self):
        '''It test that the location can check if contains another one'''
        loc_pairs = (((1, 3), (2, 5), False), ((1, 6), (4, 5), True),
                     ((2, 6), (1, 5), False), ((2, 4), (1, 5), False),
                     ((6, 10), (0, 5), False), ((6, 8), (10, 55), False))
        for loc_pair in loc_pairs:
            loc1 = Location(start=loc_pair[0][0], end=loc_pair[0][1])
            loc2 = Location(start=loc_pair[1][0], end=loc_pair[1][1])
            tuple2 = loc_pair[1]
            #it works with locations
            if loc_pair[2]:
                assert loc2 in loc1
                assert tuple2 in loc1
            else:
                #and with tuples (start, end)
                assert loc2 not in loc1
                assert tuple2 not in loc1

        #now with ints
        loc = Location(start=3, end=5)
        locs = ((1, False), (2, False), (3, True), (4, True), (5, True),
                 (6, False), (7, False))
        for int_loc in locs:
            if int_loc[1]:
                assert int_loc[0] in loc
            else:
                assert int_loc[0] not in loc

    def test_equals(self):
        '''It tests that we can know wich Locations are equal to which ones.'''
        #first by position
        loc_pairs = (((2, 5), (1, 3), False), ((2, 5), (2, 5), True),
                     ((2, 5), (2, 6), False), ((2, 5), (1, 5), False),
                     ((2, 5), (4, 5), False))
        for loc_pair in loc_pairs:
            loc1 = Location(start=loc_pair[0][0], end=loc_pair[0][1])
            loc2 = Location(start=loc_pair[1][0], end=loc_pair[1][1])
            tuple2 = loc_pair[1]
            #it works with locations
            if loc_pair[2]:
                assert loc1 == loc2
                assert loc1 == tuple2
            else:
                #and with tuples (start, end)
                assert loc1 != loc2
                assert loc1 != tuple2
        #now a problem in the other attributes
        assert Location(start=1, end=5) != \
                                   Location(start=1, end=5, strand=1)
        assert Location(start=1, end=5) != \
                                    Location(start=1, end=5, forward=False)
        assert Location(start=1, end=5) != \
                                    Location(start=1, end=5, parent='012345')
        assert Location(start=1, end=5) != \
                                    NonStaticParentLocation(start=1, end=5)

        #reverse and complement are equal
        assert Location(start=1, end=5, strand= - 1, forward=False) == \
                        Location(start=1, end=5, strand= - 1, forward=False)


class NonStaticParentLocationTests(unittest.TestCase):
    '''It tests the slicing, reverse and complement when the parent is not
    static.'''

    def test_reverse(self):
        '''Test that the Location is not mutable, when reversed it returns a
        new instance.'''
        # 012345
        #  |---|
        # |---|
        loc = NonStaticParentLocation(start=1, end=5, strand=1)
        loc2 = loc[::-1]
        assert loc2 == NonStaticParentLocation(start=0, end=4, strand=1)

        #reverse with a parent sequence of 10 residues
        seq = '0123456789'
        # 0123456789
        #  |---|
        #     |---|
        loc = NonStaticParentLocation(start=1, end=5, strand=-1, forward=False,
                                      parent=seq)
        loc2 = loc[::-1]
        assert loc2 == NonStaticParentLocation(start=4, end=8, strand=-1,
                                               forward=False)

        # 01234567890
        #      |---|
        #  |---|
        seq = '01234567890'
        loc = NonStaticParentLocation(start=5, end=9, parent=seq)
        loc2 = loc[::-1]
        assert loc2 == NonStaticParentLocation(start=1, end=5)

        try:
            #the lower than -1 strides are not supported
            loc2 = loc[::-2]
            self.fail('IndexError expected')
        except IndexError:
            pass

    def test_complement(self):
        '''Test that the Location is not mutable, when complemented it returns
        a new instance.'''
        seq = '012345'
        #with a non static parent
        loc = NonStaticParentLocation(start=1, end=5, strand=1, parent=seq)
        loc2 = loc.complement()
        assert loc2 == NonStaticParentLocation(start=1, end=5, strand=-1)

        loc = NonStaticParentLocation(start=1, end=5, strand=-1, forward=False)
        loc2 = loc.complement()
        assert loc2 == NonStaticParentLocation(start=1, end=5, strand=1,
                                               forward=False)

        seq = '01234567890'
        loc = NonStaticParentLocation(start=1, end=5, parent=seq)
        loc2 = loc.complement()
        assert loc2 == NonStaticParentLocation(start=1, end=5)

    def test_slices(self):
        '''We can get a sliced location, in a similar way as if it was a
        list.'''
        #a standard slice inside the location
        #54321 <- negative index
        #01234
        #  ---
        #  +++*
        #  012
        loc = NonStaticParentLocation(start=2, end=4, strand=1)
        loc2 = loc[2:5]
        assert loc2 == NonStaticParentLocation(start=0, end=2, strand=1)
        #same slice with a parent, the parent is not conserved
        seq = '01234'
        loc = NonStaticParentLocation(parent=seq)
        loc2 = loc[2:5]
        assert loc2 == NonStaticParentLocation(start=0, end=2)
        #with negative indexes
        loc = NonStaticParentLocation(start=2, end=4, strand=1, forward=False)
        loc2 = loc[-3:]
        assert loc2 == NonStaticParentLocation(start=0, end=2, strand=1,
                                               forward=False)
        
        #The negative indexes should work also with a parent sequence
        seq = '01234'
        loc = NonStaticParentLocation(strand=1, parent=seq)
        loc2 = loc[-3:]
        assert loc2 == NonStaticParentLocation(start=0, end=2, strand=1)

        #only partially inside the location
        #54321 <- negative index
        #01234567
        #  ---
        # +++*
        # 012
        loc = NonStaticParentLocation(start=2, end=4, strand=1)
        loc2 = loc[1:4]
        assert loc2 == NonStaticParentLocation(start=1, end=2, strand=1)
        #with negative indexes
        loc = NonStaticParentLocation(start=2, end=4, strand=1)
        loc2 = loc[-4:-1]
        assert loc2 == NonStaticParentLocation(start=1, end=2, strand=1)
    
        #01234567
        #  ---
        #   +++*
        #   01
        loc = NonStaticParentLocation(start=2, end=4, strand=1)
        loc2 = loc[3:6]
        assert loc2 == NonStaticParentLocation(start=0, end=1, strand=1)

        #01234567
        #  ---
        # +++++*
        # 0123
        loc = NonStaticParentLocation(start=2, end=4, strand=1)
        loc2 = loc[1:6]
        assert loc2 == NonStaticParentLocation(start=1, end=3, strand=1)

        #an edge case
        #01234567
        #  ---
        #    +++*
        #    0
        loc = NonStaticParentLocation(start=2, end=4, strand=1)
        loc2 = loc[4:7]
        assert loc2 == NonStaticParentLocation(start=0, end=0, strand=1)

        #an edge case
        #54321 <- negative index
        #01234
        #  ---
        #++++*
        #0123
        loc = NonStaticParentLocation(start=2, end=4, strand=1)
        loc2 = loc[:4]
        assert loc2 == NonStaticParentLocation(start=2, end=3, strand=1)
        
        #a single int as index
        #54321 <- negative index
        #01234
        #  ---
        #   ^
        #   0
        loc = NonStaticParentLocation(start=2, end=4, strand=1)
        loc2 = loc[3]
        assert loc2 == NonStaticParentLocation(start=0, end=0, strand=1)
        #with negative index
        loc2 = loc[ - 2]
        assert loc2 == NonStaticParentLocation(start=0, end=0, strand=1)

        #None returns everything
        loc = NonStaticParentLocation(start=10, end=20, strand=1)
        loc2 = loc[:None]
        assert loc2 == NonStaticParentLocation(start=10, end=20, strand=1)
        loc = NonStaticParentLocation(start=10, end=20, strand=1)
        loc2 = loc[:]
        assert loc2 == NonStaticParentLocation(start=10, end=20, strand=1)

        #a stride different than 1
        #01234567890
        #-----------
        # +   + *
        # 01234
        loc = NonStaticParentLocation(start=0, end=10, strand=1)
        try:
            loc2 = loc[1:7:4]
            self.fail('IndexError expected.')
        except IndexError:
            pass

        #negative index with a lengthier parent sequence
        #a parent sequence of 10 residues
        seq = '0123456789'
        # 0987654321 <- negative index
        # 0123456789
        #  -----
        #    -----*
        #    01234 <- new reference
        loc = NonStaticParentLocation(start=1, end=5, strand= -1, parent=seq)
        loc2 = loc[3:8]
        assert loc2 == NonStaticParentLocation(start=0, end=2, strand=-1)
        #with the negative index
        loc2 = loc[-7:-2]
        assert loc2 == NonStaticParentLocation(start=0, end=2, strand=-1)

        #index errors
        #outside the location
        loc = NonStaticParentLocation(start=10, end=20, strand=1)
        try:
            loc2 = loc[5:10]
            self.fail('IndexError expected')
        except IndexError:
            pass
        try:
            loc2 = loc[21:27]
            self.fail('IndexError expected')
        except IndexError:
            pass
        try:
            #a list would return an empty list,
            #but a biolocation is not mutable and it has not sense
            #to have an empty location
            loc2 = loc[5:5]
            self.fail('IndexError expected')
        except IndexError:
            pass
        try:
            #outside with an int
            loc2 = loc[27]
            self.fail('IndexError expected')
        except IndexError:
            pass
        try:
            #the -1 stride is not supported with start and stop
            loc2 = loc[10:20: - 1]
            self.fail('IndexError expected')
        except IndexError:
            pass
        
        #the non_static_parent does not keep the parent when sliced.
        loc = NonStaticParentLocation(parent=seq)
        loc2 = loc[2:7]
        assert loc2.parent is None

class ReparentLocationTests(unittest.TestCase):
    '''It checks the the reparent Location function behaviour'''
    def test_nonstaticparentlocation(self):
        '''It check the reparent of the NonSaticParentLocation'''
        loc = NonStaticParentLocation(start=1, end=5, strand=-1)
        seq = '0123456'
        loc2 = _reparent_location(loc, seq)
        assert loc2 == NonStaticParentLocation(start=1, end=5, strand=-1,
                                               parent=seq)

    def test_parentlocation(self):
        '''It check the reparent of the Location'''
        loc = Location(start=1, end=5, strand=-1)
        seq = '0123456'
        loc2 = _reparent_location(loc, seq)
        assert loc2 == Location(start=1, end=5, strand=-1, parent=seq)

    def test_parent_already_set(self):
        '''If the parent is already set we return the same Location.'''
        seq = '0123456'
        loc = Location(start=1, end=5, strand=-1, parent=seq)
        loc2 = _reparent_location(loc, seq)
        assert loc is loc2

class LocatableSequenceTests(unittest.TestCase):
    '''It checks the LocatableSequence behaviour'''
    def setUp(self):
        '''It creates an alignment to do some tests.'''
        #ref   0123456
        #ref   ATCTATG
        #seq1    0123
        #seq1    CTAT   Locatable
        #seq2  012345
        #seq2  tTCTag   Masked
        #seq3  43210
        #seq3  TAGAT    RevComp
        #seq4   43210
        #seq4   cGAgc   Locatable, Masked and RevComp

        # pylint: disable-msg=C0103
        self.alignment = {} #alignment seqs
        self.orig_seqs = {} #non-aligned seqs
        #parent
        ref = 'ATCTATG'
        #seq1
        seq1 = Seq('CTAT')
        self.alignment['seq1'] = locate_sequence(sequence=seq1, location=(2, 5),
                                                 strand=1, parent=ref)
        self.orig_seqs['seq1'] = seq1
        #seq2
        seq2 = Seq('tTCTag')
        self.alignment['seq2'] = locate_sequence(sequence=seq2, location=(0, 5),
                                                 strand=1, parent=ref,
                                                 mask=(1,3))
        self.orig_seqs['seq2'] = seq2
        #seq3
        seq3 = Seq('TAGAT')
        self.alignment['seq3'] = locate_sequence(sequence=seq3, location=(0, 4),
                                                 strand=-1, parent=ref,
                                                 forward=False)
        self.orig_seqs['seq3'] = seq3
        #seq4
        seq4 = Seq('cgAGc')
        def masker_x(item):
            '''It returns an x for any given item.'''
            # pylint: disable-msg=W0613
            return 'x'
        self.alignment['seq4'] = locate_sequence(sequence=seq4, location=(1, 5),
                                                 strand=-1, parent=ref,
                                                 mask=(2,3), masker=masker_x,
                                                 forward=False)
        self.orig_seqs['seq4'] = seq4

    @staticmethod
    def check_locseq(locseq, properties):
        '''It tests that the LocatableSequence has the given properties.'''
        if 'no_location' in properties:
            assert locseq.location is None
        if 'start' in properties:
            assert locseq.location.start == properties['start']
        if 'end' in properties:
            assert locseq.location.end == properties['end']
        if 'no_mask' in properties:
            assert locseq.mask is None
        if 'mask_start' in properties:
            assert locseq.mask.start == properties['mask_start']
        if 'mask_end' in properties:
            assert locseq.mask.end == properties['mask_end']
        if 'mask_parent' in properties:
            assert locseq.mask.parent is properties['mask_parent']
        if 'seq_is' in properties:
            assert locseq.sequence is properties['seq_is']
        if 'seq_equal' in properties:
            assert str(locseq.sequence) == str(properties['seq_equal'])
        if 'seq_str' in properties:
            assert str(locseq.sequence) == properties['seq_str']
        if 'strand' in properties:
            assert locseq.location.strand == properties['strand']
        if 'parent' in properties:
            assert locseq.location.parent is properties['parent']

    def test_init(self):
        '''It tests the different ways to create a LocatableSequence.'''
        #parent
        ref = 'ATCTATG'
        #seq1
        seq = 'CTAT'
        loc = NonStaticParentLocation(start=2, end=5, strand=1, parent=ref)
        locseq = LocatableSequence(sequence=seq, location=loc)
        self.check_locseq(locseq, {'start':2, 'end':5, 'no_mask':None,
                                'seq_is':seq, 'strand':1, 'parent':ref})
        #There's an easier way to do the same.
        locseq = locate_sequence(sequence=seq, location=(2, 5), strand=1,
                                  parent=ref)
        self.check_locseq(locseq, {'start':2, 'end':5, 'no_mask':None,
                                'seq_is':seq, 'strand':1, 'parent':ref})

        #seq2
        seq = 'tTCTag'
        loc = NonStaticParentLocation(start=0, end=5, strand=1, parent=ref)
        mask = NonStaticParentLocation(start=1, end=3, parent=ref)
        locseq = LocatableSequence(sequence=seq, mask=mask, location=loc)
        self.check_locseq(locseq, {'start':0, 'end':5, 'mask_start':1,
                                    'mask_end':3, 'seq_is':seq, 'strand':1,
                                    'parent':ref})
        locseq = locate_sequence(sequence=seq, location=(0, 5), strand=1,
                                  parent=ref, mask=(1,3))
        self.check_locseq(locseq, {'start':0, 'end':5, 'mask_start':1,
                                    'mask_end':3, 'seq_is':seq, 'strand':1,
                                    'parent':ref})
        
        #seq3
        seq = 'TAGAT'
        loc = NonStaticParentLocation(start=0, end=4, strand=-1, parent=ref)
        locseq = LocatableSequence(sequence=seq, location=loc)
        self.check_locseq(locseq, {'start':0, 'end':4, 'seq_is':seq,
                                    'strand':-1, 'parent':ref})
        locseq = locate_sequence(sequence=seq, location=(0, 4), strand=-1,
                                  parent=ref)
        self.check_locseq(locseq, {'start':0, 'end':4, 'seq_is':seq,
                                    'strand':-1, 'parent':ref})

        #seq4
        seq = 'cgAGc'
        mask = NonStaticParentLocation(start=2, end=3)
        def masker_x(item):
            '''It returns an x for any given item.'''
            # pylint: disable-msg=W0613
            return 'x'
        loc = NonStaticParentLocation(start=1, end=5, strand=-1, parent=ref)
        locseq = LocatableSequence(sequence=seq, location=loc, mask=mask,
                                   masker=masker_x)
        self.check_locseq(locseq, {'start':1, 'end':5, 'seq_is':seq,
                                    'strand':-1, 'parent':ref, 'mask_start':2,
                                    'mask_end':3, 'masker':masker_x})
        locseq = locate_sequence(sequence=seq, location=(1, 5), strand=-1,
                                  parent=ref, mask=(2,3), masker=masker_x)
        self.check_locseq(locseq, {'start':1, 'end':5, 'seq_is':seq,
                                    'strand':-1, 'parent':ref, 'mask_start':2,
                                    'mask_end':3, 'masker':masker_x})
        
        #a locatable sequence that only gives the start
        seq = 'cgAGc'
        locseq5 = locate_sequence(sequence=seq, location=1)
        self.check_locseq(locseq5, {'start':1, 'end':5, 'seq_is':seq,
                                    'strand':None,})
        #start at zero, the default
        locseq6 = locate_sequence(sequence=seq)
        self.check_locseq(locseq6, {'start':0, 'end':4, 'seq_is':seq})

    def test_getitem_int(self):
        '''It checks that we can get an item using an int.'''
        #seq1
        locseq1 = self.alignment['seq1']
        seq1 = self.orig_seqs['seq1']
        assert locseq1[2] == 'C'
        self.assertRaises(IndexError, locseq1.__getitem__, 1)
        self.assertRaises(IndexError, locseq1.__getitem__, 6)
        assert len(locseq1) == len(seq1)
        #seq2
        locseq2 = self.alignment['seq2']
        seq2 = self.orig_seqs['seq2']
        assert locseq2[1] == 'T'
        self.assertRaises(IndexError, locseq2.__getitem__, 0)
        self.assertRaises(IndexError, locseq2.__getitem__, 5)
        self.assertRaises(IndexError, locseq2.__getitem__, 6)
        assert len(locseq2) == len(seq2)
        #seq3
        locseq3 = self.alignment['seq3']
        seq3 = self.orig_seqs['seq3']
        assert locseq3[0] == 'A'
        assert locseq3[4] == 'A'
        self.assertRaises(IndexError, locseq3.__getitem__, 5)
        assert len(locseq3) == len(seq3)
        #seq4
        locseq4 = self.alignment['seq4']
        seq4 = self.orig_seqs['seq4']
        assert locseq4[2] == 'C'
        assert locseq4[1] == 'x'
        self.assertRaises(IndexError, locseq4.__getitem__, 0)
        assert locseq4[5] == 'x'
        self.assertRaises(IndexError, locseq4.__getitem__, 6)
        assert len(locseq4) == len(seq4)

    def test_getitem_slice(self):
        '''We test that we get the correct LocatableSequence after a slice'''
        #now we slice it
        #ref   0123456
        #ref   ATCTATG
        #seq1    0123
        #seq1    CTAT   Locatable
        #seq2  012345
        #seq2  tTCTag   Masked
        #seq3  43210
        #seq3  TAGAT    RevComp
        #seq4   43210
        #seq4   cGAgc   Locatable, Masked and RevComp
        #      0123456
        #slice  [----[
        #new    012345
        #ref    TCTATG
        #seq1    0123
        #seq1    CTAT   Locatable
        #seq2   01234
        #seq2   TCTag   Masked
        #seq3   3210
        #seq3   AGAT    RevComp
        #seq4   43210
        #seq4   cGAgc   Locatable, Masked and RevComp
        locseq1 = self.alignment['seq1']
        seq1 = self.orig_seqs['seq1']
        locseq11 = locseq1[1:6]
        self.check_locseq(locseq11, {'start':1, 'end':4, 'no_mask':None,
                                'seq_is':seq1, 'strand':1})
        locseq2 = self.alignment['seq2']
        locseq21 = locseq2[1:6]
        self.check_locseq(locseq21, {'start':0, 'end':4, 'mask_start':0,
                                'mask_end':2, 'seq_equal':'TCTag'})
        locseq3 = self.alignment['seq3']
        locseq31 = locseq3[1:6]
        self.check_locseq(locseq31, {'start':0, 'end':3, 'no_mask':None,
                                'seq_equal':'TAGA', 'strand':-1})
        locseq4 = self.alignment['seq4']
        seq4 = self.orig_seqs['seq4']
        locseq41 = locseq4[1:6]
        self.check_locseq(locseq41, {'start':0, 'end':4, 'mask_start':2,
                                'mask_end':3, 'seq_is':seq4, 'strand':-1})

        #no support for slices different than 1 or -1
        self.assertRaises(IndexError, locseq1.__getitem__, slice(None, None, 3))

    def test_reverse(self):
        '''It tests that we can reverse the LocatableSequence.'''
        #reverse
        #ref   0123456
        #ref   ATCTATG
        #seq1    0123
        #seq1    CTAT   Locatable
        #seq2  012345
        #seq2  tTCTag   Masked
        #seq3  43210
        #seq3  TAGAT    RevComp
        #seq4   43210
        #seq4   cGAgc   Locatable, Masked and RevComp
        #reversed
        #ref   0123456
        #seq1   0123
        #seq1   TATC
        #seq2   012345
        #seq2   gaTCTt
        #seq3    43210
        #seq3    TAGAT
        #seq4   43210
        #seq4   cgAGc
        locseq1 = self.alignment['seq1']
        locseq1r = locseq1[::-1]
        self.check_locseq(locseq1r, {'start':1, 'end':4, 'no_mask':None,
                                'seq_equal':'TATC', 'strand':1})
        locseq2 = self.alignment['seq2']
        locseq2r = locseq2[::-1]
        self.check_locseq(locseq2r, {'start':1, 'end':6, 'mask_start':2,
                               'mask_stop':4, 'seq_equal':'gaTCTt', 'strand':1})
        locseq3 = self.alignment['seq3']
        locseq3r = locseq3[::-1]
        self.check_locseq(locseq3r, {'start':2, 'end':6, 'no_mask':None,
                                'seq_equal':'TAGAT', 'strand': - 1})
        locseq4 = self.alignment['seq4']
        locseq4r = locseq4[::-1]
        self.check_locseq(locseq4r, {'start':1, 'end':5, 'mask_start':1,
                               'mask_stop':2, 'seq_equal':'cGAgc', 'strand':-1})

    def test_complement(self):
        '''It tests that we can complement a LocatableSequence.'''
        #complement
        #ref   0123456
        #ref   ATCTATG
        #seq1    0123
        #seq1    CTAT   Locatable
        #seq2  012345
        #seq2  tTCTag   Masked
        #seq3  43210
        #seq3  TAGAT    RevComp
        #seq4   43210
        #seq4   cGAgc   Locatable, Masked and RevComp
        #complemented
        #ref   0123456
        #ref   ATCTATG
        #ref   TAGATAC
        #seq1    0123
        #seq1    GATA 
        #seq2  012345
        #seq2  tTCTag   Masked
        #seq2  aAGAtc
        #seq3  43210
        #seq3  ATCTA
        #seq4   43210
        #seq4   gCTcg
        locseq1 = self.alignment['seq1']
        locseq1c = locseq1.complement()
        self.check_locseq(locseq1c, {'start':2, 'end':5, 'no_mask':None,
                                     'seq_equal':'GATA', 'strand':1})
        locseq2 = self.alignment['seq2']
        locseq2c = locseq2.complement()
        self.check_locseq(locseq2c, {'start':0, 'end':5, 'mask_start':1,
                                     'mask_stop':3, 'seq_equal':'aAGAtc'})
        locseq3 = self.alignment['seq3']
        locseq3c = locseq3.complement()
        self.check_locseq(locseq3c, {'start':0, 'end':4, 'no_mask':None,
                                     'seq_equal':'ATCTA', 'strand': - 1})
        locseq4 = self.alignment['seq4']
        locseq4c = locseq4.complement()
        self.check_locseq(locseq4c, {'start':1, 'end':5, 'mask_start':2,
                                     'mask_stop':3, 'seq_equal':'gcTCg',
                                     'strand':-1})
        #The sequence should know how to be complemented
        locseqerr = locate_sequence(sequence='ACTGTC')
        self.failUnlessRaises(AttributeError, locseqerr.complement)

    def test_as_maskedsequence(self):
        '''It tests it's behaviour as a MaskedSequence.'''
        #seq  0123
        #mask x--x
        seq = [0, 1, 2, 3]
        mask = locate_sequence(sequence=seq, mask=(1, 2))
        self.check_locseq(mask, {'mask_start':1, 'mask_stop':2, 'seq_is':seq,
                                  'mask_parent':seq})

        #a mask is not given
        #     0123456
        #seq  no_mask
        #mask -------
        seq = 'no_mask'
        mask = locate_sequence(sequence=seq)
        self.check_locseq(mask, {'no_mask':None, 'seq_is':seq})

        #     01234
        #seq  ACTGT
        #mask x--xx
        #with a masker function
        seq = 'ACTGT'
        def lower(item):
            '''It returns a lower str'''
            return item.lower()
        mask = locate_sequence(sequence=seq, mask=(1, 2), masker=lower)
        self.check_locseq(mask, {'mask_start':1, 'mask_stop':2, 'seq_is':seq,
                                  'mask_parent':seq})

        #check mask outside the sequence
        #     012345
        #seq  ACTGT
        #mask x-----x
        seq = 'ACTGT'
        try:
            locate_sequence(sequence=seq, mask=(1, 5))
            self.fail('We expected a ValueError')
        except ValueError:
            pass

        #slicing int
        #seq  0123
        #mask x--x
        seq = [0, 1, 2, 3]
        mask = locate_sequence(sequence=seq, mask=(1, 2))
        assert mask[1] == 1
        self.assertRaises(IndexError, mask.__getitem__, 0)

        #masker
        #     01234
        #seq  ACTGT
        #mask x--xx
        #with a masker function
        seq = 'ACTGT'
        mask = locate_sequence(sequence=seq, mask=(1, 2),
                               masker=lambda x: x.lower())
        assert mask[0] == 'a'
        assert mask[1] == 'C'

        #slicing slice
        #masker
        #     0123456
        #seq  ACTGTGT
        #mask x----xx
        #with a masker function
        seq = 'ACTGTGT'
        mask = locate_sequence(sequence=seq, mask=(1, 4))
        #      0123456
        #seq   ACTGTGT
        #mask  x----xx
        #slice x---
        mask2 = mask[:4]
        assert str(mask2.sequence) == str(seq[:4])
        assert mask2.mask.start == 1
        assert mask2.mask.end == 3
        self.check_locseq(mask2, {'mask_start':1, 'mask_stop':3,
                                  'seq_equal':seq[:4]})
        #reverse
        #      0123456
        #seq   ACTGTGT
        #mask  x----xx
        #seq   TGTGTCA
        #mask  xx----x
        mask3 = mask[::-1]
        self.check_locseq(mask3, {'mask_start':2, 'mask_stop':5,
                                  'seq_equal':seq[::-1]})

        #      0123456789
        #seq   ACTGTGT
        #mask  x----xx
        #slice         --
        seq = 'ACTGTGT'
        mask = locate_sequence(sequence=seq, mask=(1, 4))
        self.assertRaises(IndexError, mask.__getitem__, slice(8, 10))

        #slices different than 1 or -1 are not supported.
        self.assertRaises(IndexError, mask.__getitem__, slice(0, 5, 3))

    def test_as_revcomp(self):
        '''It checks the behaviour as a reverse and complementary container.'''
        #slicing int
        
        #        0123
        #rev     TGAC
        #revcomp ACTG
        #mutable sequence reverse
        revseq = locate_sequence(Seqmut('CAGT'), forward=False)
        assert revseq[0] == 'T'
        assert revseq[3] == 'C'
        #mutable sequence reverse and complementary
        revseq = locate_sequence(Seqmut('CAGT'), strand=-1, forward=False)
        #it will raise a ValueError because we don't support complementing
        #of mutable sequences. We would be changing a sequence outside the
        #LocatableSequence
        self.assertRaises(ValueError, revseq.__getitem__, 0)
        #non mutable sequence reverse
        revseq = locate_sequence(Seq('CAGT'), forward=False)
        assert revseq[0] == 'T'
        assert revseq[3] == 'C'
        #non mutable sequence reverse and complementary
        revseq = locate_sequence(Seq('CAGT'), strand=-1, forward=False)
        assert revseq[0] == 'A'
        assert revseq[3] == 'G'
        #if a sequence non-complementable we get an AttributeErrror
        revseq = locate_sequence('wrong', strand=-1)
        self.assertRaises(AttributeError, revseq.__getitem__, 0)
        
        #slicing slice
        #        0123
        #revcomp ACTG
        #slice   --
        #non mutable sequence
        revseq = locate_sequence(Seq('CAGT'), strand=-1, forward=False)
        revseq2 = revseq[:2]
        assert str(revseq2.sequence) == 'GT'
        #mutable sequence
        revseq = locate_sequence(Seqmut('CAGT'), strand=-1, forward=False)
        revseq2 = revseq[:2]
        assert str(revseq2.sequence) == 'GT'
        #a plain sequence can't be complemented
        #        0123
        #revcomp ACTG
        #rev     TGAC
        #    seq GTCA
        #        3210
        #slice   ---------
        revseq = locate_sequence('CAGT', strand=-1, forward=False)
        revseq2 = revseq[:8]
        assert str(revseq2.sequence) == 'CAGT'

        #with negative indexes
        #plain sequence
        #       -4321
        #        0123
        # revseq CAGT
        #    seq GTCA
        #        3210
        #slice    --
        revseq = locate_sequence(Seqmut('CAGT'), strand=-1, forward=False)
        revseq2 = revseq[-3:-1]
        assert str(revseq2.sequence) == 'AG'

        #it raises an IndexError outside the boundaries
        revseq = locate_sequence(Seqmut('CAGT'), strand=-1, forward=False)
        self.assertRaises(IndexError, revseq.__getitem__, 4)
        self.assertRaises(IndexError, revseq.__getitem__, -5)
        self.assertRaises(IndexError, revseq.__getitem__, slice(4, 5))
        self.assertRaises(IndexError, revseq.__getitem__, slice(-6, -5))

    @staticmethod
    def test_outside_limits_config():
        '''It checks that we can configure what we get outside the seq limits'''
        import biolib.contig as contig
        locseq = contig.locate_sequence('ACTG', location=7, mask=(2, 3))
        #we can configure what we get outside the sequence
        contig.empty_region_seq_builder = lambda : ''
        assert locseq[0] == ''
        #now what we get in the masked region
        contig.default_masker = lambda x: 'x'
        assert locseq[8] == 'x'

    @staticmethod
    def test_str():
        'It test the str writting of the LocatableSequence.'
        #with an start
        locseq = locate_sequence(sequence='ACTG', location=4)
        assert locseq.__str__() == '    ACTG'
        #with a parent
        locseq = locate_sequence(sequence='ACTG', location=4,
                                 parent='01234567890')
        assert locseq.__str__() == '    ACTG   '
        #with a mask and a masker
        locseq = locate_sequence(sequence='ACTG', mask=(1, 2),
                                 masker=lambda x: x.lower())
        assert locseq.__str__() == 'aCTg'
        #with a mask and no masker
        locseq = locate_sequence(sequence='ACTG', mask=(1, 2))
        assert locseq.__str__() == 'xCTx'
        #with everything
        locseq = locate_sequence(sequence='ACTG', location=1, mask=(1, 2),
                                 masker=lambda x: x.lower(), parent='012345')
        assert locseq.__str__() == ' aCTg '
        
class ContigTests(unittest.TestCase):
    '''It checks the Contig class behaviour.'''
    def test_basic_behaviour(self):
        '''It tests the init and the basic slicing.'''
        #a basic assembly with all sequences starting at the same place
        #aligment
        #     0123456
        #seq1 ATCTGCA
        #seq2 ATGTGCA
        #seq3 ATGTGCT
        seq = 'ATCTGCA'
        alignment = Contig(sequences=[seq, 'ATGTGCA', 'ATGTGCT'])
        assert len(alignment) == 3
        assert alignment.ncols == 7
        #some slicing
        assert alignment[0, 0] == 'A'
        assert alignment[0] is seq
        assert alignment[1, 1:3] == 'TG'
        #now we get a column
        assert alignment[1:, 6] == 'AT'

        #the same with seqRecords
        seqrec1 = SeqRecord('ATCTGCA')
        seqrec2 = SeqRecord('ATGTGCA')
        alignment = Contig(sequences=[seqrec1, seqrec2, SeqRecord('ATGTGCT')])
        assert len(alignment) == 3
        assert alignment.ncols == 7
        #if we want the SeqRecod.seq we should ask for it
        assert alignment[0, 0].seq == 'A'
        #a seqrecord slice returns a seqrecord slice
        assert alignment[0] is seqrec1
        assert alignment[1, 1:3].seq == 'TG'
        #now we get a column
        assert alignment[:, 2].seq == 'CGG'

        #now a contig without consensus
        #ref   012345
        #seq1    0123
        #seq1    CTAT   Locatable
        #seq2  012345
        #seq2  tTCTag   Masked
        #seq3  43210
        #seq3  TAGAT    RevComp
        #seq4   43210
        #seq4   cGAgc   Locatable, Masked and RevComp

        contig = Contig()
        #seq0
        seq = 'TAGAT'
        contig.append(seq)
        #seq1
        seq = 'CTAT'
        contig.append_to_location(seq, start=2)
        #seq2
        seq = 'tTCTag'
        mask = (1, 3)
        contig.append_to_location(seq, mask=mask)
        #seq3
        seq = 'TAGAT'
        contig.append_to_location(seq, strand=-1)
        #seq4
        seq = 'cgAGc'
        mask = (2, 3)
        contig.append_to_location(seq, start=1, strand=-1, mask=mask)

        assert len(contig) == 5
        assert contig.ncols == 6
        assert contig[1, 2] == 'C'
        self.assertRaises(IndexError, contig.__getitem__, (1, 0))
        #if we have located a sequence inside the contig with an slice we'll
        #get a LocatableSequences
        assert contig[0, 1:3].sequence == 'AG'
        assert contig[1, 1:3].sequence == 'C'

        #now a contig composed by seqwithquality
        #ref   012345
        #seq1    ATCT fwd
        #qual1   5678
        #seq2  GCTA   revcomp
        #qual2 89123
        contig = Contig()
        seq1 = SeqWithQuality(Seq('ATCT'), [5, 6, 7, 8])
        contig.append_to_location(seq1, start=2)
        seq2 = SeqWithQuality(Seq('ATCG'), [3, 2, 1, 9])
        contig.append_to_location(seq2, strand=-1, forward=False)

        assert len(contig) == 2
        assert contig.ncols == 6
        assert contig[1, 0].seq == 'C'
        assert contig[0, 1:4].sequence.seq == 'AT'
        assert contig[0, 1:4].sequence.qual[0] == 5

        #now we get a column
        col = contig[:, 2]
        #the seqwithqual
        assert col.seq == 'AA'
        assert str(col.qual) == str([5, 2])
        #only the seq
        assert contig[:, 0].seq == 'C'
        #only the qual
        assert str(contig[:, 5].qual) == str([8])

        #let's slice an assembly
        #      0123456789
        #seq1  ATCT
        #seq2     TGACA
        #seq3         ACT
        contig2 = Contig()
        contig2.append('ATCT')
        contig2.append_to_location('TGACA', start=3)
        contig2.append_to_location('ACT', start=7)
        contig21 = contig2[:, 1:7]
        assert len(contig21) == 2
        assert contig21.ncols == 6
        assert contig21[0, 0] == 'T'
        assert contig21[1, 2] == 'T'

    def test_consensus(self):
        '''It tests that a consensus can be added.'''
        contig1 = Contig(consensus='ACGT')
        assert contig1.consensus == 'ACGT'

        #now we slice the contig
        contig2 = contig1[:, :2]
        assert contig2.consensus == 'AC'

    def test_empty_seqs_slice(self):
        '''We want also the empty seqs, could we get them?'''
        seq = [locate_sequence('ACTG', location=5), 'AAAA', 'TTTTTT']
        contig = Contig(seq)
        # Slice slice combination
        #default behaviour
        new_contig = contig[:, 0:1]
        assert len(new_contig) == 2
        #now we want the empty seqs
        contig.return_empty_seq = True
        new_contig = contig[:, 0:1]
        assert new_contig[0] is None
    
        # int slice combination
        contig.return_empty_seq = True
        new_contig = contig[0, 2:3]
        assert new_contig is None
        
    def test_empty_seqs_int(self):
        '''We want also the empty seqs, could we get them?'''
        seq = [locate_sequence('ACTG', location=5), 'AAAA', 'TTTTTT']
        contig = Contig(seq)        

        # These are with slice int combination
        #default behaviour
        contig.return_empty_seq = False
        new_contig = contig[:, 1]
        assert len(new_contig) == 2
        #now we want the empty seqs
        contig.return_empty_seq = True
        new_contig = contig[:, 0]
        assert new_contig[0] is None
        
        #These are with int int combination
        #default behaviour
        contig.return_empty_seq = False
        try:
            #pylint: disable-msg=W0104
            #the statement does has an effect
            contig[0, 1]
            self.fail("Index Error expected")
        except IndexError:
            pass
        
    @staticmethod
    def test_consensus_empty():
        ''' We test here if the contig watches what returns the contig when 
        you ask a column that doesn't exist'''
        consensus = locate_sequence('TT', location=0)
        contig = Contig(sequences=['AAAA'], consensus = consensus)
        assert  contig[:, 4:4].consensus  == None

    @staticmethod
    def test_contig_str():
        'It checks the str contig representation'
        consensus = locate_sequence('AA', location=1)
        contig = Contig(sequences=['xAAx'], consensus = consensus)
        assert str(contig) == '0              ->xAAx\nconsensus      -> AA\n'
   
def main():
    '''It runs the tests'''
    unittest.main()

if __name__ == '__main__':
    main()
