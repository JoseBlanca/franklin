'''This module provides the code to represent a sequence contig.'''
from biolib.locatable_sequence import locate_sequence, slice_to_range

#for configuration I like lowercase variables
#pylint: disable-msg=C0103
def default_masker_function(seq):
    '''It returns an index error to be used as default masker'''
    #pylint: disable-msg=W0613
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

    def __str__(self):
        'It produces a string that shows all the reads and consensus'
        tostring = []
        #all the reads
        for index, read in enumerate(self):
            #maybe the sequence has a name property
            try:
                name = read.name
            except AttributeError:
                #if it's a LocatableSequence the name should be the in
                #the sequence property
                try:
                    name = read.sequence.name
                except AttributeError:
                    name = str(index)
            tostring.append(name.ljust(20))
            tostring.append('->')
            tostring.append(str(read))
            tostring.append('\n')
        #the consensus
        if self.consensus is not None:
            tostring.append('consensus'.ljust(20))
            tostring.append('->')
            tostring.append(str(self.consensus))
            tostring.append('\n')
        return ''.join(tostring)

    def __repr__(self):
        '''It writes  the reads of the contigs'''
        info_repr = ''
        for read in self:
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
    
    def _getitem_slice_slice(self, row_index, col_index):
        '''It returns a new Contig with the sliced rows.'''
        #we create the new assembly that will be returned
        new_assembly = self.__class__()
        for row_int in slice_to_range(row_index, len(self)):
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
        for row_int in slice_to_range(row_index, len(self)):
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
