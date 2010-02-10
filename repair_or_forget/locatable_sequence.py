'''It provides a LocatableSequence and a Location class

LocatableSequence is a wrapper around a sequence that allows to mask and
change the coordinate system of a sequence without changing the sequence
itself.

Location represents a segment on a sequence.
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

from biolib.biolib_utils import get_start_end

def raise_error_outside_limits():
    '''It returns an index error to be used outside the seqs boundaries'''
    raise IndexError('Index outside sequence limits')
#by default when an item is requested outside the sequence boundaries an
#IndexError is raised
empty_region_seq_builder = raise_error_outside_limits

def locate_sequence(sequence, location=None, mask=None, masker=None,
                    strand=None, forward=True, parent=None):
    '''It creates a new LocatableSequence.
    
    It creates a LocatableSequence with the given location and mask.
    keyword arguments:
    sequence -- The sequence to be located
    location -- an tuple of ints with the position of the first and the
                last bases in the parent sequence or an int with just the
                start (default None)
    mask -- A tuple with the first an last non-masked bases (default None)
    masker -- The masker function (default None)
    strand -- The location strand (default None)
    forward -- The location direction (default True)
    parent -- The parent sequence to which the LocatableSequence will be 
              refered (default None)
    '''
    #strs to ints
    if mask:
        mask = [int(number) for number in mask]
    # pylint: disable-msg=R0913
    #I know that it has a lot of arguments, but the LocatableSequence is
    #a complex class
    loc = None
    if location is None:
        start = 0
        end   = len(sequence) - 1
    elif isinstance(location, int):
        start = location
        end   = location + len(sequence) - 1
    else:
        start = int(location[0])
        end   = int (location[1])
    loc = NonStaticParentLocation(start=start, end=end, strand=strand,
                                  forward=forward, parent=parent)
    mask_instance = None
    if mask is not None:
        mask_instance = NonStaticParentLocation(start=mask[0], end=mask[1],
                                                parent=sequence)
    return LocatableSequence(sequence=sequence, location=loc,
                             mask=mask_instance, masker=masker)

def _sequence_end(sequence):
    '''If there is a parent where does the parent ends, else
    self.end.'''
    #if the sequence has an stop, that's what we need.
    #We don't want just the length, because the sequence could
    #start at a position different than zero   
    try:
        return sequence.end
    except AttributeError:
        return len(sequence) - 1

def _positive_int(index, parent, sequence):
    '''It returns the same int index, but positive.'''
    if index is None:
        return None
    elif index < 0:
        if parent is not None:
            length = _sequence_end(parent)
        else:
            length = _sequence_end(sequence)
        return length + index + 1
    return index
def _positive_slice(index, parent, sequence):
    '''It returns the same slice but with positive indexes.'''
    return slice(_positive_int(index.start, parent, sequence),
                 _positive_int(index.stop, parent, sequence),
                 index.step)

def _positive(index, sequence, parent=None):
    '''Given an index it returns its positive version.
    
    The index coordinate system can be referenced to a sequence or to its parent.
    If the parent is given it will be used to calculate the positive index over
    the sequence.
    
            0123456789 <- positive
            9876543210 <- negative index
    parent  aparentseq
    sequece aseq
    index    876 <- negative
             123 <- positive
    The end of the parent and the sequence will be calculated using the end 
    property, if it the sequence doesn't have an end property the length will
    be used instead. The reason for this is that not all sequences start at zero.
    keyword arguments:
    index -- the index to transform (an int or an slice)
    sequence -- the sequence in which the index is
    parent -- the parent in which the index is (default None)
    '''
    if isinstance(index, slice):
        return _positive_slice(index, parent, sequence)
    else:
        return _positive_int(index, parent, sequence)
def _move_index(index, amount):
    '''It returns a new index displaced an amount of units.
    
    The index can be None, int or slice.
    '''
    if index is None:
        return None
    elif isinstance(index, int):
        return index + amount
    else:
        start = index.start
        stop = index.stop
        if start is None:
            new_start = None
        else:
            new_start = start + amount
        if stop is None:
            new_stop = None
        else:
            new_stop = stop + amount
        return slice(new_start, new_stop, index.step)

class _SequenceAttribute(object):
    '''This class is used to create the properties that return different
    properties of the contained sequence (self.seq).
    It's used to provided the syntactic sugar: locseq.seq[1] locseq.qual[1]
    '''
    #This is just syntax sugar, it should do very little
    # pylint: disable-msg=R0903
    def __init__(self, parent, req_attr):
        '''The initialitation.
        keyword arguments:
        parent -- The parent object that owns the property (e.g. the locseq)
        req_attr -- The property requested to the contained sequence.
        '''
        self._parent = parent
        self._req_attr = req_attr
    def __getitem__(self, index):
        '''It returns an item or an slice of the requested property.'''
        #_sequence_getitem is private, but this is not standard client code,
        #it's almost internal to LocatableSequence
        #pylint: disable-msg=W0212
        return self._parent._sequence_getitem(index=index,
                                             attr=self._req_attr)

class LocatableSequence(object):
    '''It's a container for sequence like objects.

    It accepts lists, strings, Seqs, etc.
    It places the given sequence in relation with a parent sequence. It's
    main aim is to locate a read inside a contig. Besides changing the
    coordinate system it supports adding a mask at the beginning and at the
    end of the sequence.
    The objects derived from this class are non-mutable.
    '''
    #it has noly one public method, but a complex pythonic behaviour
    # pylint: disable-msg=R0903
    def __init__(self, sequence, location=None, mask=None, masker=None):
        '''To initialize an object only a sequence-like object is required.

        If a Location is not given a default one will be created with the
        start at 0, end at length-1 and strand 0.
        If no mask is given, the whole sequence will be unmasked.
        A masker function could be given. This function would transform the
        items from the masked region. When this function is given no IndexError
        will be raised in the masked regions. The function should take an
        item like the ones given by the sequence.
        keyword arguments:
        sequence -- The sequence-like object stored (e.g. Seq, list, str)
        location -- A Location with the coords in the ref seq (default None).
        mask -- A Location with the coords in the seq of the masked region.
        masker -- A function that returns a masked item.
        '''
        self._sequence = sequence
        self._location = None
        self._set_location(location)
        self._mask = None
        self._set_mask(mask)
        self._masker = masker
        #the special seq and qual properties
        self.seq  = _SequenceAttribute(self, 'seq')
        self.qual = _SequenceAttribute(self, 'qual')
        #some caches
        self._masked_comp_seq_dict_cache = None
        self._masked_comp_seq_cache = None

    @staticmethod
    def _check_location(seq, location):
        '''It checks that the sequence fits into the space provided by the
        location.
    
        If it not fits it will raise an error.
        '''
        if location is None:
            return
        if  len(seq) != len(location):
            msg = 'Sequence and Location lengths do not match.'
            raise ValueError(msg)

    def _set_location(self, location):
        '''It sets the location and it checks that the seqs fits in it.
        
        If location is None it will create one.
        '''
        if location is None:
            start = 0
            end = len(self.sequence) - 1
            location = NonStaticParentLocation(start=start, end=end)
        else:
            self._check_location(self.sequence, location)
        self._location = location
 
    def _set_mask(self, mask):
        '''It sets the mask and it checks that it fits in the sequence.'''
        # check_location_fits
        if mask is not None and\
           not mask in Location(start=0, end=len(self.sequence) - 1):
            raise ValueError('The mask does not fit in the sequence')
        if mask is not None:
            mask = _reparent_location(location=mask, parent=self.sequence)
        self._mask = mask

    def _get_sequence(self):
        '''It returns the sequence instance.'''
        return self._sequence
    sequence = property(_get_sequence)

    def _get_location(self):
        '''It returns the location.'''
        return self._location
    location = property(_get_location)

    def _get_mask(self):
        '''It returns the mask.'''
        return self._mask
    mask = property(_get_mask)

    def _get_masker(self):
        '''It returns the masker.'''
        return self._masker
    masker = property(_get_masker)

    def _masked_comp_seq_dict(self):
        '''It complements, reverses and masks self.sequence.
        
        It stores the result in a dict with three dicts that correspond to the first
        masked region, the non-masked sequence and the second masked region.
        Every dict has three keys start, end and sequence. sequence can be
        None in the case of the masked regions.
        '''
        if self._masked_comp_seq_dict_cache is not None:
            return self._masked_comp_seq_dict_cache
        #this fucntion does all the transformations to the sequence but the
        #change in the coord system.
        location = self.location
        strand   = location.strand
        forward  = location.forward
        mask     = self.mask
        seq      = self.sequence
        end      = len(seq) - 1
        #what do we need to calculate?
        #  mask_1_region | sequence | mask_2_region
        #  start     end  start  end  start     end
        m1_start = None
        m1_end   = None
        m1_seq   = None
        seq_start = 0
        seq_end  = end
        seq      = seq
        m2_start = None
        m2_end   = None
        m2_seq   = None
        #complement
        if strand == -1:
            #The complementary sequence
            type_before_complement = type(seq)
            seq = seq.complement()
            #maybe the sequence was mutable, in that case the sequence now
            #will be None. We raise an error in this case because we don't 
            #want to mess with sequences that could be being used outside
            #this class for other porpouses
            if type(seq) != type_before_complement:
                msg = 'The sequence seems mutable, we support only non-mutables'
                raise ValueError(msg)
        #mask
        if mask is not None:
            m1_start  = 0
            seq_start = mask.start
            m1_end    = seq_start - 1
            seq_end   = mask.end
            m2_start  = seq_end + 1
            m2_end    = len(seq) - 1
            m1_seq    = seq[:seq_start]
            m2_seq    = seq[m2_start:]
            seq       = seq[seq_start:m2_start]
        #reverse
        if not forward:
            #now the mask1 is mask2 and mask2 is mask1,
            #and the indexes should be recalculated counting from the end
            if m2_end:
                length = m2_end + 1
            else:
                length = seq_end + 1
            if m1_start is not None:
                m2_start  = length - 1 - m1_start
                m2_end    = length - 1 - m1_end
            if m2_start is not None:
                m1_start  = length - 1 - m2_start
                m1_end    = length - 1 - m2_end
            seq_start = length - 1 - seq_start
            seq_end   = length - 1 - seq_end
            seq       = seq[::-1]
            old_m2_seq = m2_seq
            if m1_seq:
                m2_seq = m1_seq[::-1]
            if m2_seq:
                m1_seq = old_m2_seq[::-1]
        #now we can return
        res = {}
        res['mask1'] = {}
        res['seq']   = {}
        res['mask2'] = {}
        res['mask1']['start'] = m1_start
        res['mask1']['end']   = m1_end
        res['mask1']['seq']   = m1_seq
        res['seq']['start']   = seq_start
        res['seq']['end']     = seq_end
        res['seq']['seq']     = seq
        res['mask2']['start'] = m2_start
        res['mask2']['end']   = m2_end
        res['mask2']['seq']   = m2_seq
        self._masked_comp_seq_dict_cache = res
        return res

    def _check_in_ok_region(self, index, error_in_out, error_in_mask):
        '''It raises an error if we're in an inconvenient region.
        
        The index is in the self coord system
        '''
        seq = self._masked_comp_seq_dict()
        location = self.location
        start = location.start
        end   = location.end
        if isinstance(index, slice):
            index_start = index.start
            index_end   = index.stop - 1
        else:
            index_start = index
            index_end   = index
        #we check if we're in a region that should raise an error
        #are we outside the seq?
        if (error_in_out and
            (index_start, index_end) not in Location(start, end)):
            raise IndexError('Sequence asked outside the region covered by seq')
        #are we in the mask
        if error_in_mask:
            index_loc = Location(index_start, index_end)
            #where does the umasked region starts and ends?
            if seq['mask1']['start'] is not None:
                mask_start = seq['mask1']['start'] + start
                mask_end   = seq['mask1']['end']   + start
                if index_loc.overlaps((mask_start, mask_end)):
                    msg = 'The requested region covers the mask'
                    raise IndexError(msg)
            #where does the umasked region starts and ends?
            if seq['mask2']['start'] is not None:
                mask_start = seq['mask2']['start'] + start
                mask_end   = seq['mask2']['end']   + start
                if index_loc.overlaps((mask_start, mask_end)):
                    msg = 'The requested region covers the mask'
                    raise IndexError(msg)

    def _masked_comp_seq(self):
        '''It returns self.sequence masked and complemented'''
        if self._masked_comp_seq_cache is not None:
            return self._masked_comp_seq_cache
        seq = self._masked_comp_seq_dict()

        joined = None
        #the mask
        masker = self.masker
        if seq['mask1']['start'] is not None:
            if masker is None:
                #the sequence unmasked
                joined = seq['mask1']['seq'][:]
            else:
                joined = masker(seq['mask1']['seq'])
        #the sequence
        if joined is None:
            joined = seq['seq']['seq']
        else:
            joined += seq['seq']['seq']
        #the mask2
        if seq['mask2']['start'] is not None:
            if masker is None:
                joined += seq['mask2']['seq']
            else:
                joined += masker(seq['mask2']['seq'])
        self._masked_comp_seq_cache = joined
        return joined

    def _sequence_getitem(self, index, attr=None):
        '''A get item for int indexes function that takes into account the
        mask and the need to complement.
        
        It requires an index in the coord system of self.
        It takes into account if we need the complement or not.
        The difference with __getitem__ is that this one returns an instance
        like self.sequence, not a LocatableSequence
        '''
        #we take into account the space from 0 till the seq start
        pseudo_seq_index = _move_index(index, -self.location.start)
        seq = self._masked_comp_seq_dict()

        if self.mask and not self.masker:
            error_in_mask = True   #error if we're in the mask
        else:
            error_in_mask = False
        error_in_out  = True    #error if we're outside the seq
        self._check_in_ok_region(index, error_in_out, error_in_mask)

        #we join the dict in one piece
        seq =  self._masked_comp_seq()
        #do we have to return directly from seq or from one of its attributes?
        if attr is None:
            seq_attr = seq
        else:
            seq_attr = getattr(seq, attr)
        return seq_attr[pseudo_seq_index]

    def __str__(self):
        'It returns the LocatableSequence as a string'
        #we ask for everything from o to the end of the parent
        seq = self._masked_comp_seq_dict()

        #now we get the complete sequence in str format
        location = self.location
        start    = location.start
        parent   = location.parent

        tostring = []
        #the LocatableSequence has several sections
        #first part till the start
        if start > 0:
            tostring.append(' ' * start)
        #the mask
        masker = self.masker
        if seq['mask1']['end'] is not None:
            if masker is None:
                tostring.append(' ' * len(seq['mask1']['seq']))
            else:
                tostring.append(str(masker(seq['mask1']['seq'])))
        #the sequence
        tostring.append(str(seq['seq']['seq']))
        #the mask2
        if seq['mask2']['start'] is not None:
            if masker is None:
                tostring.append(' ' * len(seq['mask2']['seq']))
            else:
                tostring.append(str(masker(seq['mask2']['seq'])))
        #the last part from sequence end till the parent end
        if parent is not None:
            if seq['mask2']['end'] is not None:
                seq_end = seq['mask2']['end']
            else:
                seq_end = seq['seq']['end']
            last_empty_spaces = len(parent) - seq_end - start - 1
            if last_empty_spaces:
                tostring.append(' ' * last_empty_spaces)
        seq_str = ''.join(tostring)

        return seq_str

    def __repr__(self):
        ''' It writes'''
        return self._location.__repr__() + self._sequence.__repr__() + "\n"
        
    def __getitem__(self, index):
        '''It returns an item or a new LocatableSequence.

        The index should be in the coord. system of self.
        '''
        location = self.location
        #we need the index in the coord system of self and of self.seq
        self_index = index
        #we have to be sure that it is positive
        index = _positive(index, self, location.parent)
        #we transform the index in the coord. system of self.seq
        #if self.seq is reversed the index should be moved one position to the
        #left, because otherwise it won't represent the same region. stop is 
        #different than end by one, that's the problem
        #self coord                0123456
        #region covered by index     xxxx
        #direct index                [   [ start=2, stop=6
        #reversed index             ]   ]  stop=1, start=5 
        if not location.forward and isinstance(index, slice):
            seq_index = location.get_location_index(_move_index(index, -1))
        else:
            seq_index = location.get_location_index(index)
        if seq_index is None:
            #we're outside the sequence, by default we raise an IndexError
            return empty_region_seq_builder()
        #do we return an item?
        if isinstance(index, int):
            return self._sequence_getitem(index)
        #or we return a new LocatableSeq?
        else:
            #we do not support step slices different than None, 1 or -1
            if self_index.step not in (None, 1, -1):
                raise IndexError('slice step not supported (only 1, -1, None)')
            #we slice everything
            new_loc = None
            if location is not None:
                new_loc = location[self_index]
            mask = self.mask
            new_mask = None
            if mask is not None:
                new_mask = mask[seq_index]
            #if the slice covers the whole sequence and is not a reverse
            #we don't need to copy (slice) the sequence
            start = seq_index.start
            stop = seq_index.stop
            #we get the requested seq property, self.sequence,
            #self.sequence.qual or whatever
            sequence = self.sequence
            if (start is None or start == 0) and \
               (stop is None or stop >= len(sequence)) and \
               seq_index.step != -1:
                new_seq = sequence
            else:
                new_seq = sequence[seq_index]
            return self.__class__(sequence=new_seq, location=new_loc,
                                              mask=new_mask, masker=self.masker)

    def __len__(self):
        '''It returns the length.'''
        return len(self.sequence)

    def complement(self):
        '''It returns a new LocatableSequence with a complemented seq.

        This is not reverse complement, so the original seq and the new seq
        wont't be the same.
        '''
        #the sequence should know how to complement itself
        new_seq = self.sequence.complement()
        new_loc = self.location
        if new_loc is not None:
            new_loc = new_loc[:]    #we need a copy becase Location is inmutable
        new_mask = self.mask
        if new_mask is not None:
            new_mask = new_mask[:]
        return self.__class__(sequence=new_seq, location=new_loc, mask=new_mask,
                              masker=self.masker)

    empty_region_seq_builder = raise_error_outside_limits

    def __getattr__(self, attrname):
        '''It returns the attrs of the located sequence.'''
        #there are two kinds of attributes, the ones that should be returned
        #straight from self.sequence like name and the ones that should take
        #into account the different indexing between self and self.sequence and
        #the reverse, complement and masking like seq and qual
        #those attrs that have to change the index are not handle by this
        #method
        return getattr(self.sequence, attrname)

    def add_mask(self, mask, masker=None):
        '''It returns a new masked LocatableSequence.
        
        The new mask applied to the sequence will be the intersection beteween
        thre previous mask and the new one.
        If a masker is given it will replace the one present in the Location.
        '''
        loc = self.location
        if loc.start != 0 or not loc.forward or loc.strand == -1:
            #for that to work the given mask should be converted to the
            #coord system of the sequence
            raise NotImplementedError('Remasking with location not ready yet')
        if masker is None:
            masker = self.masker
        old_mask = self.mask
        if old_mask is not None:
            new_mask = old_mask.intersection(mask)
        else:
            new_mask = mask
        return self.__class__(self.sequence, mask=new_mask, masker=masker,
                             location=loc)


def _reparent_location(location, parent):
    '''Given a Location it returns a new Location with the given parent.
    
    This function is required because Location is not mutable and it is a
    common thing to ask for a new Location with the same start, end, strand,
    and forward, but with a different parent.
    If the parent in the given Location is the same as the given parent it
    will return the given location not a new one. 
    '''
    #has the Location already the parent?
    if location.parent is parent:
        return location
    return location.__class__(start=location.start, end=location.end,
                              strand=location.strand, forward=location.forward,
                              parent=parent)

class Location(object):
    '''A Location represents a sequence segment.

    It represents a fragment of a sequence.
    Its basic use is to relate its own coordinate system with the
    sequence one.
         012345678901234678901234567890
         thisisjustanexampleofasequence
                0123456789012
                <-alocation->
    Optionally it can be related to a parent sequence. That would mean that 
    the  Location is a fragment of that parent sequence.
    The Location can be sliced, reversed and complemented. Since is a
    non-mutable class any of these operations would return  a new Location.
    If it has a parent it assumes that this parent is static.
    An static parent will not suffer the slicing, reverse and complement, only
    the Location would. This is important to calculate the relation between
    the Location and the new parent.
    The only requirement to be a Location parent is to have a length. 
    '''
    def __init__(self, start=None, end=None, strand=None, forward=True,
                 parent=None):
        '''It initializes the Location.

        The Location can be initialized with an start and end coordinates or
        with a parent. In the latter case if the end is not given il will be
        set to len(parent) - 1.
        keyword arguments:
        start -- the begining of the Location (default 0).
        end -- the end of the Location (default len(parent) - 1)
        strand -- The strand can be 1 or -1. (default None)
        forward -- The Location can be in the same direction as the parent or
                   reversed (default True).
        parent -- The parent sequence (default None)
        static_parent -- would the parent not be affected by the slicing?
                         (default True)
        '''
        #According to pylint this class should has only 5 arguments, but
        #I think that, in this case, they are justified
        # pylint: disable-msg=R0913
        self._parent = None
        self._set_parent(parent)
        self._strand = None
        self._set_strand(strand)
        self._start = None
        self._end = None
        self._set_start_end(start, end)
        self._forward = None
        self._set_forward(forward)

    def _set_parent(self, parent):
        '''It sets the parent.'''
        self._parent = parent
    def _get_parent(self):
        '''It returns the parent.'''
        return self._parent
    parent = property(_get_parent)

    def _set_start_end(self, start, end):
        '''It sets the start and end of the Location.'''
        parent = self.parent
        if start is None:
            self._start = 0
        else:
            self._start = start
        if end is None:
            if parent is None:
                raise ValueError('end cannot be None if parent is None.')
            else:
                self._end = len(parent) - 1
        else:
            self._end = end
        #some checks
        if self.end < self.start:
            raise ValueError ('start should be lower than end.')
        if self.start < 0:
            raise ValueError ('start cannot be negative.')

    def _get_start(self):
        '''It returns the Location start.'''
        return self._start
    start = property(_get_start)
    def _get_end(self):
        '''It returns the Location end.'''
        return self._end
    end = property(_get_end)

    def _set_forward(self, forward):
        '''It sets wether the Location is in the same direction as the
        parent or not.'''
        self._forward = forward
    def _get_forward(self):
        '''It returns wether the Location is in the same direction as the
        parent or not.'''
        return self._forward
    forward = property(_get_forward)

    def _set_strand(self, strand):
        '''It sets the strand Location.

        It should be 1, -1 or None.
        '''
        if strand not in (1, None, -1):
            raise ValueError('Acceptable strand values are: 1, None and -1.')
        self._strand = strand
    def _get_strand(self):
        '''It returns the strand Location.'''
        return self._strand
    strand = property(_get_strand)

    def __len__(self):
        '''It returns the Location len.

        Mind that the Location can start at a position different than 0.
        '''
        return self.end - self.start + 1
    def __repr__(self):
        ''' It prints the location start and end'''
        return "Location: start  %d, end %d\t" % (self._start, self._end)
    def apply_to_parent(self):
        '''It returns the parent sequence fragment that corresponds to the
        Location.

        For instance:
        >>> parent_seq = 'ACATGCTA'
        >>> loc = Location(start=1, end=4, parent=parent_seq)
        >>> print loc.apply_to_parent()
        CATG
        '''
        if self.parent is None:
            raise RuntimeError('The parent is not defined.')
        stop = self.end + 1
        answer = self.parent[slice(self.start, stop)]
        if not self.forward:
            answer = answer[::-1]
        if self.strand == -1:
            answer = answer.complement()
        return answer

    def get_parent_index(self, index):
        '''Given an index in the Location coordinate system it returns the
        same index in the parent coordinate system.

        The index can be an int or an slice.
        '''
        def parent_index_int(index):
            '''Given an index in the Location coordinate system it returns
            the same index in the parent coordinate system.

            The index should be an int.
            '''
            if index is None:
                return None
            if self.forward:
                return index + self.start
            else:
                return self.end - index
        if index is None:
            return None
        if isinstance(index, int):
            return parent_index_int(index)
        #the index is a slice
        return slice(parent_index_int(index.start),
                     parent_index_int(index.stop), index.step)

    def get_location_index(self, index):
        '''Given an index in the parent coordinate system it returns the
        same index in the Location coordinate system.

        The index can be an int or an slice.
        '''
        self_forward = self.forward
        self_start = self.start
        self_end = self.end

        def loc_index_int(index):
            '''Given an index in the parent coordinate system it returns the
            same index in the Location coordinate system.

            The index should be an int.
            '''
            if index is None:
                return None
            if self_forward:
                new_index = index - self.start
            else:
                new_index = self_end - index
            #if we're outside the Location it should be None
            if new_index < 0 or new_index > (self_end - self_start):
                new_index = None
            return new_index
        if index is None:
            return None
        if isinstance(index, int):
            return loc_index_int(index)
        #new index
        new_start = loc_index_int(index.start)
        new_stop = loc_index_int(index.stop)
        #the index is a slice
        if self_forward:
            return slice(new_start, new_stop, index.step)
        else:
            return slice(new_stop, new_start, index.step)

    def _reversed_location(self):
        '''It retusns a revesed Location.'''
        return self.__class__(start=self.start, end=self.end,
                    strand=self.strand, forward=not(self.forward),
                    parent=self.parent)

    def _sliced_location(self, start, end, o_start=None):
        '''It returns an sliced Location.
        
        It uses the slice start and end and optionally the original slice
        start.
        '''
        #I'm aware of the unused argument, but this arguments is required
        #by NonStaticParentLocation
        # pylint: disable-msg=W0613

        if self.forward:
            return self.__class__(start=start, end=end, strand=self.strand,
                                  forward=True, parent=self.parent)
        else:
            end2 = self.end - (start - self.start)
            start2 = self.start - (end - self.end)
            return self.__class__(start=start2, end=end2, strand=self.strand,
                                  forward=False, parent=self.parent)

    def __getitem__(self, index):
        '''It returns a new sliced Location.

        The index should be in the coordinate system of the parent sequence.
        The parent sequence is not sliced, only the Location is sliced.
                       before         after
           seq_coords  0123456        0123456
           parent_seq  ACTGCAT        ACTGCAT
           Loc_coords    0123           01
           Location      xxxx   ->      xx
           Slice        [--[
        The parent of the new Location is the same as the old one.
        '''
        #positive
        index = _positive(index, self, self.parent)
        #the int index is treated as a slice of length 1
        if isinstance(index, int):
            s_start = index
            s_stop = index + 1
            sr_step = 1
        else:
            s_start = index.start
            s_stop = index.stop
            sr_step = index.step
        #step should be 1, None or - 1
        if sr_step not in (1, None, -1):
            raise IndexError('Only 1, None, -1 steps are allowed.')
        #slice and reverse is not allowed at the same time
        if sr_step == -1 and (s_start is not None or s_stop is not None):
            raise IndexError('slicing an reverse not allowed at the same time.')
        #do we want to reverse
        if sr_step == -1:
            return self._reversed_location()

        #the slice should overlap with the location
        if s_stop is None:
            s_end = None
        else:
            s_end = s_stop - 1
        if not self.overlaps((s_start, s_end)):
            raise IndexError('Index not covered by the Location')
        #both start and end should be inside the location to continue with the
        #calculations
        if s_start is None or self.start >= s_start:
            sr_start = self.start
        else:
            sr_start = s_start
        if s_end is None or self.end <= s_end:
            sr_end = self.end
        else:
            sr_end = s_end

        #different ways of calculate the results if parent is static or not
        return self._sliced_location(sr_start, sr_end, s_start)

    def complement(self):
        '''It returns a new complemented location.

        It just changes the strand.
        If the parent is not-static the new Location wouldn't know who its
        parent is.
        '''
        new_strand = self.strand
        if new_strand:  #different than 0
            new_strand = -new_strand
        return self.__class__(start=self.start, end=self.end,
                              strand=new_strand, forward=self.forward,
                              parent=self.parent)

#    @staticmethod
#    def _get_start_end(location):
#        '''It accepts an int, Location or tuple and it returns the start, end,
#        forward and strand.'''
#        #int
#        if isinstance(location, int):
#            start = location
#            end = location
#        #tuple
#        elif isinstance(location, tuple):
#            start = location[0]
#            end = location[1]
#        #location
#        else:
#            start = location.start
#            end = location.end
#        return start, end

    def overlaps(self, location):
        '''It returns True if the locations overlap, False otherwise.

        It can accepts Locations, tuples with the start and end or int.
        Keyword argument:
        location -- A Location to compare.
        '''
        start, end = get_start_end(location)
        #self    -----     -----     -----     --------
        #range  ---      ---------    ------    -----
        #if start or end is None we convert them to ints to allow the
        #comparisons
        #this is an optimization based on a cprofiler run
        self_start = self.start
        self_end = self.end
        if start is None:
            start = 0
        if end is None:
            end = self_end
            #if start was bigger than self.end and end was None
            #at this point end could be lower than start
            if end < start:
                end = start
        if self_start >= start and self_start <= end:
            return True
        if self_end >= start and self_end <= end:
            return True
        if self_start <= start and self_end >= end:
            return True
        return False

    def __contains__(self, location):
        '''It returns True if the given location is totally inside self.

        It can accepts Locations, tuples with the start and end or int.
        Keyword argument:
        location -- A Location or an int to compare.
        '''
        start, end = get_start_end(location)
        #slice
        #self  ------------
        #range   --------- 
        if start >= self.start and end <= self.end:
            return True
        return False

    @staticmethod
    def parent_is_static():
        '''The parent of this class is static, so we return True.'''
        return True

    def __eq__(self, location):
        '''It returns True if this location is equal to the given location.

        It can accepts Locations, tuples with the start and end or int.
        Keyword argument:
        location -- A Location to compare.
        '''
        #self  ---------
        #range --------- 
        start, end = get_start_end(location)
        #here we also need strand, forward and static_parent  
        try:
            strand = location.strand
        except AttributeError:
            strand = None
        try:
            forward = location.forward
        except AttributeError:
            forward = True
        try:
            parent = location.parent
        except AttributeError:
            parent = None
        try:
            static = location.parent_is_static()
        except AttributeError:
            static = None

        answer = True
        if self.start != start:
            answer = False
        if self.end != end:
            answer = False
        if self.strand != strand:
            answer = False
        if self.parent != parent:
            answer = False
        if self.forward != forward:
            answer = False
        if static is not None and static != self.parent_is_static():
            answer = False
        return answer
    def intersection(self, location2):
        ''' It returns a new location with the intersection of the two locations
         given'''
        forward  = self.forward
        strand   = self.strand
        forward2 = location2.forward
        strand2  = location2.strand

        if forward != forward2 or strand != strand2 :
            raise ValueError("Strand and fordward must be the same for both")
            
        start1 = self.start
        end1  = self.end
        start2 = location2.start
        end2  = location2.end
        
        
        if start1 < start2 and end1 < start2:
            return None
        elif start1 > end2 and end1 > end2:
            return None
        else:
            if start1 > start2:
                start = start1
            else:
                start = start2
            if end1 > end2:
                end = end2
            else:
                end = end1
        
        return self.__class__(start=start, end=end, strand=strand, \
                              forward=forward, parent=self.parent)
         
class NonStaticParentLocation(Location):
    '''Its a Location with a parent affected by the slicing.
 
    The parent will be also sliced, complemented and reversed.
    This is important to calculate the relation between the Location and the
    new parent.
    The only requirement to be a Location parent is to have a length. 
    '''
    def __init__(self, start=None, end=None, strand=None, forward=True,
                 parent=None):
        '''It initializes the Location.

        The Location can be initialized with an start and end coordinates or
        with a parent. In the latter case if the end is not given il will be
        set to len(parent) - 1.
        keyword arguments:
        start -- the begining of the Location (default 0).
        end -- the end of the Location (default len(parent) - 1)
        strand -- The strand can be 1 or -1. (default None)
        forward -- The Location can be in the same direction as the parent or
                   reversed (default True).
        parent -- The parent sequence (default None)
        static_parent -- would the parent not be affected by the slicing?
                         (default True)
        '''
        #According to pylint this class should has only 5 arguments, but
        #I think that, in this case, they are justified
        # pylint: disable-msg=R0913

        Location.__init__(self, start=start, end=end, strand=strand,
                          forward=forward, parent=parent)

    @staticmethod
    def parent_is_static():
        '''The parent of this class is static, so we return True.'''
        return True

    def _reversed_location(self):
        '''It retusns a revesed Location.'''
        if self.parent is not None:
            r_length = len(self.parent)
        else:
            r_length = self.end + 1
        start = r_length - self.end - 1
        end = r_length - self.start - 1
        #the new forward won't change because the parent is also
        #being reversed. The parent is not conserved because we don't
        #now who the new parent will be (it will be the result of the
        #reverse of the old parent
        return self.__class__(start=start, end=end, strand=self.strand,
                              forward=self.forward)

    def _sliced_location(self, start, end, o_start=None):
        '''It returns an sliced Location'''
        if o_start is None:
            o_start = 0
        offset = start - o_start
        start2 = offset
        end2 = end - start + offset
        return self.__class__(start=start2, end=end2, strand=self.strand,
                              forward=self.forward)

    def complement(self):
        '''It returns a new complemented location.

        It just changes the strand.
        If the parent is not-static the new Location wouldn't know who its
        parent is.
        '''
        new_strand = self.strand
        if new_strand:  #different than 0
            new_strand = -new_strand
        return self.__class__(start=self.start, end=self.end,
                              strand=new_strand, forward=self.forward)

def slice_to_range(index, length):
    '''Given a slice and the length of the sequence it returns a range.'''
    try:
        start = index.start
    except AttributeError:
        return [index]
    if start is None:
        start = 0
    stop = index.stop
    if stop is None:
        stop = length
    step = index.step
    if step is None:
        step = 1
    return range(start, stop, step)
