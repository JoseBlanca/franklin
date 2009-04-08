'''This module provides the code to represent a sequence contig.'''

#for configuration I like lowercase variables
#pylint: disable-msg=C0103
def default_masker_function(seq):
    '''It returns an index error to be used as default masker'''
    raise IndexError('Index in masked sequence region')
#by default when a sequence in a masked region is requested an IndexError is
#raised
default_masker = default_masker_function

def raise_error_outside_limits():
    '''It returns an index error to be used outside the seqs boundaries'''
    raise IndexError('Index outside sequence limits')
#by default when an item is requested outside the sequence boundaries an
#IndexError is raised
empty_region_seq_builder = raise_error_outside_limits

class Contig(object):
    '''It represents a list of aligned sequences with a list-like interface.

    It can be used for clustal-like alignments or for contigs. For the contigs
    the sequences should be introduced as LocatableSequence.
    The sequences can be str, Seq, Qual, list or whatever sequence-like objects.
    The biopython SeqRecord is not yet a sequence-like because it lacks the
    __getitem__ method, but this class also supports it because it knows that
    it should access its seq property.
    '''
    def __init__(self, sequences=None, consensus=None):
        '''Initialize a new Contig.

        It accepts a list or a generator of sequence-like objects. If not given
        these sequences could be appended afterwards.
        keyword arguments:
        sequences -- a list or iterator with sequence like objects (default=[])
        consensus -- a consensus sequence (default None)
        '''
        self._seqs = []
        self._consensus = None
        #the Contig is located within the start and end
        #if it's a clustal-like assembly start would be 0 and end len -1
        self._start = 0
        self._end = 0
        #we have to control that the attrs should be only set once.
        self._attrs_set = False
        #By default we asume that we're dealing with an alignment not with a
        #contig
        self._contig_mode = False
        #by default the contig won't return empty lines, but that behavior
        #can be modified. If return_empty_seq is True it will return rows
        #with a None value when a sequence slicing returns an IndexError
        self.return_empty_seq = False

        if sequences is not None:
            self.extend(sequences)
        self.consensus = consensus

    def _set_consensus(self, consensus):
        '''It sets the consensus sequence.'''
        self._consensus = consensus
    def _get_consensus(self):
        '''It returns the consensus sequence.'''
        return self._consensus
    consensus = property(_get_consensus, _set_consensus)

    def extend(self, sequences):
        '''It extends the list by appending elements from the iterable.'''
        for seq in sequences:
            self.append(seq)

    def append(self, sequence):
        '''It appends a sequence to the end of the sequence list.'''
        self.append_to_location(sequence)

    def append_to_location(self, sequence, start=0, strand=1, forward=True,
                           mask=None, masker=None):
        '''It appends a sequence to the end of the sequence list.

        The sequence added can start at 0 (the default) or at any other position
        specified by start. E.g. in a contig not all sequences start at the
        same position.
        assem 0123456
        seq1  ATCTA   start = 0
        seq2   TCTATT start = 1
        Also the sequence can be in the reverse-complemented strand if strand is
        -1. Acceptable values for strand are 1, or -1.
        assem 0123456
        seq1  ATCTA   start = 0
        seq2   TCTATT start = 1 strand = -1 seq = 'AATAGA'
        A fragment at the begining and at the end can be masked. The region not
        masked should be specified by a tuple with two ints (start and end)
        using the coordinate system of the read. The start should be the first
        non masked letter and the end the last one. By default everything is
        unmasked.
        The default behaviour is to raise an IndexError when a masked position
        is requested, but this can be changed with a masker function. The
        masker function should take an item (nucleotide or quality) and return
        the masked value.
        assem 01234567890
        seq1     aATCTAt  mask=(1, 5) sequence='AATCTAT'
                          masker=lambda item: item.lower()
        Be aware that when start, strand or mask are used the sequence will be
        wrapped in a LocatableSequence and that is what it will be returned
        by the Contig.__getitem__ method.
        '''
        # pylint: disable-msg=R0913
        #This is a helper function and it really requieres those optional
        #arguments, I don't know how to remove them
        
        # pylint: disable-msg=W0511
        #TODO. check alphabets
        #if the sequences have an alphabet check that the alphabets are
        #compatible.
        def wrap_seqs():
            '''It wraps all the existing sequences into LocatableSeqs.'''
            for i, seq in enumerate(self._seqs):
                seq = locate_sequence(sequence=seq, parent=self)
                self._seqs[i] = seq
        #are we dealing with a contig or with an alignment?
        if not self._contig_mode and \
           (start != 0 or strand == -1 or not forward or mask is not None):
            #if there were already sequences added to the contig we have to
            #wrap them into LocatableSeqs
            wrap_seqs()
            self._contig_mode = True
        #if we're in contig mode we wrap the sequences into LocatableSequences
        if self._contig_mode:
            sequence = locate_sequence(sequence=sequence, location=start,
                                       strand=strand, forward=forward,
                                       parent=self, mask=mask, masker=masker)
        #we update the start and end of the alignment
        #where does this sequence starts and ends?
        try:
            loc = sequence.location
            start = loc.start
            end = loc.end
        except AttributeError:
            start = 0
            if sequence is not None:
                end = len(sequence) - 1
            else:
                end = 0
        if start < self._start:
            self._start = start
        if end > self._end:
            self._end = end
        #we add the sequence to the list
        self._seqs.append(sequence)

    def __len__(self):
        '''It returns the number of rows in the alignment.'''
        return len(self._seqs)

    def __get_ncols(self):
        '''It returns the number of columns in the alignment.'''
        return self._end - self._start + 1
    ncols = property(__get_ncols)

    def __repr__(self):
        ''' It writes  the reads of the contigs'''
        info_repr = ''
        for read in self._seqs:
            info_repr += read.__repr__()
        return "Consesus:" + self._consensus.__repr__() + info_repr
        
    def _getitem_int_slice(self, row_index, col_index):
        '''It returns a row or an sliced row, (a sliced seq  or a seq).'''
        row = self._seqs[row_index]
        #do we really want a slice or we want the whole row?
        if col_index.start is None and col_index.stop is None and \
           col_index.step is None:
            return row
        else:
            try:
                return row[col_index]
            except IndexError:
                #some sequences might not span in this column range
                if self.return_empty_seq:
                    return  None
                else:
                    raise                 
    @staticmethod            
    def _slice_to_range(index, length):
        '''Given a slice and the length of the sequence it returns a range.'''
        start = index.start
        if start is None:
            start = 0
        stop = index.stop
        if stop is None:
            stop = length
        step = index.step
        if step is None:
            step = 1
        return range(start, stop, step)

    def _getitem_slice_slice(self, row_index, col_index):
        '''It returns a new Contig with the sliced rows.'''
        #we create the new assembly that will be returned
        new_assembly = self.__class__()
        for row_int in self._slice_to_range(row_index, len(self)):
            # pylint: disable-msg=W0704
            #It's ok to do nothing here
            try:
                row = self._getitem_int_slice(row_int, col_index)
                new_assembly.append(row)
            except IndexError:
                #some sequences might not span in this column range
                if self.return_empty_seq:
                    new_assembly.append(None)
                else:
                    pass
        #now we have to slice the consensus
        try:
            consensus = self.consensus
            if consensus is not None:
                new_assembly.consensus = consensus[col_index]
        except IndexError:
            if self.return_empty_seq:
                new_assembly.consensus = None
            else:
                pass
                
        return new_assembly
    def _getitem_slice_int(self, row_index, col_index):
        '''It returns a column.'''
        #which row range are we requesting?
        items    = None
        seq_like = None
        for row_int in self._slice_to_range(row_index, len(self)):
            try:
                item = self._seqs[row_int][col_index]
            except IndexError:
                #this row does not cover this column
                if self.return_empty_seq:
                    item = None
                else:
                    continue
            #are we dealing with seq-like items like strs or Seq
            #or with non-seq items like the int of a quality
            if seq_like is None:
                try:
                    len(item)
                    seq_like = True
                except TypeError:
                    seq_like = False
            #if is a seq_like we sum the items
            if seq_like:
                if items is None:
                    items = item
                else:
                    items += item
            else:
            #else, we create a list
                if items is None:
                    items = [item]
                else:
                    items.append(item)
        return items

    def __getitem__(self, index):
        '''It returns a sequence or a new Contig.

        The index can be an item or a tuple of two items, these two items
        correspond to the row and the column indexes. These items can be 
        int or slice.
        Depending on the index requested different items will be returned.
        int   -> a row (a sequence) like the one that was appended.
        slice -> a new Contig with the sequences contained in the slice.
        None, int -> a column wrapped in the class of the kind that form the
                     rows.
        None, slice -> a new Contig with all the row sequences sliced.
        int, int -> An item in the assembly.
        int, slice -> a sequence slice (the same as assembly[int][slice])
        slice, int -> A column with not all sequences taken into account.
        slice, slice -> A new Contig with the sequences contained in the slice
                        sliced.
        '''

        #fisrt we make sure that we have a tuple as index
        if type(index) in (int, slice):
            row_index = index
            col_index = None
        else:
            row_index = index[0]
            col_index = index[1]
        #now we make sure that none of them is None
        if row_index is None:
            row_index = slice(None, None, None)
        if col_index is None:
            col_index = slice(None, None, None)
        #at this point there's only four options (int, int), (int, slice),
        #(slice, int) or (slice, slice)
        row_type = type(row_index)
        col_type = type(col_index)
        if row_type == int and col_type == int:
            try:
                return self._seqs[row_index][col_index]
            except IndexError:
                if self.return_empty_seq:
                    return None
                else:
                    raise
        elif row_type == int and col_type == slice:
            return self._getitem_int_slice(row_index, col_index)
        elif row_type == slice and col_type == int:
            return self._getitem_slice_int(row_index, col_index)
        elif row_type == slice and col_type == slice:
            return self._getitem_slice_slice(row_index, col_index)


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
    end_cache = None
    def sequence_end(sequence, end_cache=end_cache):
        '''If there is a parent where does the parent ends, else
        self.end.'''
        #if the sequence has an stop, that's what we need.
        #We don't want just the length, because the sequence could
        #start at a position different than zero
        if end_cache is not None:
            return end_cache
        try:
            end_cache = sequence.end
        except AttributeError:
            end_cache = len(sequence) - 1
        return end_cache
    def positive_int(index):
        '''It returns the same int index, but positive.'''
        if index is None:
            return None
        elif index < 0:
            if parent is not None:
                length = sequence_end(parent)
            else:
                length = sequence_end(sequence)
            return length + index + 1
        return index
    def positive_slice(index):
        '''It returns the same slice but with positive indexes.'''
        return slice(positive_int(index.start), positive_int(index.stop),
                     index.step)

    if isinstance(index, slice):
        return positive_slice(index)
    else:
        return positive_int(index)
def _move_index(index, amount):
    '''It returns a new index displaced an amount of units.
    
    The index can be None, int or slice.
    '''
    if index is None:
        return None
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

    def _sequence_property(self, complement=False):
        '''It returns the self.sequence or its complement.'''
        # pylint: disable-msg=W0201
        #is defined here, but is used only here
        self._comp_cache = {}
        def complement_cache(sequence):
            '''It returns a complemented sequence.
            
            It keeps a cache of complemented sequences to speed up the int
            slicing. It has no sense to complement the whole sequence every
            time we ask for an item in the complemented sequence
            '''
            instance_id = id(sequence)
            if not instance_id in self._comp_cache:
                self._comp_cache[instance_id] = sequence.complement()
            return self._comp_cache[instance_id]
        sequence = self.sequence
        #do we want the complementary?
        if complement:
            type_before_complement = type(sequence)
            sequence = complement_cache(sequence)
            #maybe the sequence was mutable, in that case the sequence now
            #will be None. We raise an error in this case because we don't 
            #want to mess with sequences that could be being used outside
            #this class for other porpouses
            if type(sequence) != type_before_complement:
                msg = 'The sequence seems mutable, we support only non-mutables'
                raise ValueError(msg)
        #Which property do we return? self.sequence or self.sequence.seq?
        return sequence
        
    def _masked_comp_seq_int(self, index):
        '''A get item for int indexes function that takes into account the
        mask and the need to complement.
        It requires an int index.
        It takes into account if we need the complement or not.
        '''
        #This is the sequence that we should return, it can be
        #self.sequence, self.sequence.seq, self.sequence.qual or whatever
        complement = False
        if self.location.strand == -1:
            complement = True
        sequence = self._sequence_property(complement)
        #if there is no mask we just return using the seq getitem
        if self._mask is None: 
            return sequence.__getitem__(index)
        #index is int
        if index in self.mask:
            return sequence.__getitem__(index)
        else:
            if self.masker is not None:
                return self.masker(sequence.__getitem__(index))
            else:
                return default_masker(sequence.__getitem__(index))

    def __repr__(self):
        ''' It writes'''
        return self._location.__repr__() + self._sequence.__repr__() + "\n"
        
    def __getitem__(self, index):
        '''It returns an item or a new LocatableSequence.

        The index should be in the coord. system of self.
        '''
        #we need the index in the coord system of self and of self.seq
        self_index = index
        #we have to be sure that it is positive
        index = _positive(index, self, self.location.parent)
        #we transform the index in the coord. system of self.seq
        #if self.seq is reversed the index should be moved one position to the
        #left, because otherwise it won't represent the same region. stop is 
        #different than end by one, that's the problem
        #self coord                0123456
        #region covered by index     xxxx
        #direct index                [   [ start=2, stop=6
        #reversed index             ]   ]  stop=1, start=5 
        if not self.location.forward and isinstance(index, slice):
            seq_index = self._location.get_location_index(
                                                  _move_index(index, -1))
        else:
            seq_index = self._location.get_location_index(index)
        if seq_index is None:
            #we're outside the sequence, by default we raise an IndexError
            return empty_region_seq_builder()
        #do we return an item?
        if isinstance(seq_index, int):
            return self._masked_comp_seq_int(seq_index)
        #or we return a new LocatableSeq?
        else:
            #we do not support step slices different than None, 1 or -1
            if self_index.step not in (None, 1, -1):
                raise IndexError('slice step not supported (only 1, -1, None)')
            #we slice everything
            loc = self.location
            new_loc = None
            if loc is not None:
                new_loc = loc[self_index]
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

    @staticmethod
    def _get_start_end(location):
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

    def overlaps(self, location):
        '''It returns True if the locations overlap, False otherwise.

        It can accepts Locations, tuples with the start and end or int.
        Keyword argument:
        location -- A Location to compare.
        '''
        start, end = self._get_start_end(location)
        #self    -----     -----     -----     --------
        #range  ---      ---------    ------    -----
        #if start or end is None we convert them to ints to allow the
        #comparisons
        if start is None:
            start = 0
        if end is None:
            end = self.end
        if self.start >= start and self.start <= end:
            return True
        if self.end >= start and self.end <= end:
            return True
        if self.start <= start and self.end >= end:
            return True
        return False

    def __contains__(self, location):
        '''It returns True if the given location is totally inside self.

        It can accepts Locations, tuples with the start and end or int.
        Keyword argument:
        location -- A Location or an int to compare.
        '''
        start, end = self._get_start_end(location)
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
        start, end = self._get_start_end(location)
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
