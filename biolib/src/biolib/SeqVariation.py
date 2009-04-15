'''
Created on 2009 mar 25

@author: peio
'''
from biolib.contig import NonStaticParentLocation, Contig

class _SeqVarConf(object):
    '''This class contains some switches to configure to your needs
    '''
    #No public methods at all, only some properties.
    #we could use a namedtuple but it would be less clear.
    #pylint: disable-msg=R0903
    def __init__(self, min_num_of_reads=2, only_snp=False,
                indel_char='-', empty_char=''):
        ''' Here we initialize teh object with the configuration we want.
        
        This configuration will affect to all classes and functions from this
        method.
        min_num_of_reads is the number of times that an alleles has to be read
        to be considered. This is used to ease the removal of some sequence
        errors.
        If only_snp is set to True the indels will be ignored.
        keyword arguments:
            min_num_of_reads: an int (default 2)
            only_snp: a bool (default False)
            indel_char: the character used for the inserts (default '-')
        '''
        self.min_num_of_reads = min_num_of_reads
        self.only_snp         = only_snp
        self.indel_char       = indel_char
        self.empty_char       = empty_char


CONFIG = _SeqVarConf()

class SeqVariation(object):
    '''
    This class is used to represent any kind of sequence variation in 
    and alignment. The variation can be a Snp or a InDel, or ...
    '''

    def __init__(self, alleles, name=None, location=None, alignment=None):
        '''
        This class collect all the information about each secuence variation
        in an alignment. Sequence variations could be SNPs or InDels. There is a
        filter in the module configuration to use only SNPs in this class.
        To initialize the only required information is the alleles.
        Keyword arguments:
            alleles - . The aleles contains a dict with the alleles and the
                        number of times that it has been read.
                        example:
                            alleles= {'A':1, 'T':4}
            name     - Name of the variation (default None)
            location - Location of the SeqVariation. Location class object or
                       int (default None).
            alignemt - Alignment where it procedes. It can be a Contig class 
                       object (default None).
        '''
        self.name         = name
        self._num_reads   = None  # the alleles dict with the alleles as keys
        self._set_alleles(alleles)
        self.location     = location
        self.alignment    = alignment
    
    def __repr__(self):
        ''' It prints a readable SeqVariation information '''
        
        pstring  = 'SequenceVariation: ' + str(self.name) + '\n'
        if self.location:
            pstring += 'Location:' + self.location.__repr__() + '\n'
#        if self.alignment:
#            pstring += self.alignment.name + '\n'
        pstring += self._num_reads.__repr__()
        return pstring
    
    def _get_alleles(self):
        ''' It returns a dict with the aleles.'''
        return self._num_reads
    
    def _set_alleles(self, alleles):
        ''' It sets the alleles dict.
        It takes into account the min_number_of_reads and it remove the ones
        that do not meet that criteria.
        '''
        for allele in alleles.keys():
            #we remove the alleles with not enough reads
            if isinstance(alleles[allele], int):
                allele_count = alleles[allele]
            else:
                allele_count = len(alleles[allele])
            if allele_count < CONFIG.min_num_of_reads:
                del alleles[allele]
                continue
            #we remove the indels if we don't want them
            if CONFIG.only_snp and allele == CONFIG.indel_char:
                del alleles[allele]
        self._num_reads = alleles
    alleles = property(_get_alleles)

    def is_complex(self):
        '''It returns True if the variation is a complex '''
        if self.kind() == 'complex':
            return True
        else:
            return False
        
    def is_indel(self):
        '''It returns True if the variation is an indel.'''
        if self.kind() == 'indel':
            return True
        else:
            return False
        
    def is_snp(self):
        '''It returns True if the variation is a snp '''
        if self.kind() == 'snp':
            return True
        else:
            return False
    
    def kind(self):
        '''It returns the kind of variation: snp, indel or complex.'''
        num_alleles    = len(self._num_reads)
        inchar = CONFIG.indel_char
        #if there is only one allele is invariable 
        if num_alleles < 2:
            var_kind = 'invariable'
        elif num_alleles == 2:
            allele1, allele2 = self._num_reads.keys()
            if allele1 == inchar * len(allele1) and inchar not in allele2:
                var_kind = 'indel'
            elif allele2 == inchar * len(allele2) and inchar not in allele1:
                var_kind =  'indel'
            elif inchar in allele1 or inchar in allele2:
                var_kind =  'complex'
            else:
                var_kind =  'snp'
        elif num_alleles >2:
            alleles = self._num_reads.keys()
            if len(alleles[0]) == 1 and inchar not in alleles:
                var_kind = 'snp'
            else:
                var_kind =  'complex'
        return var_kind

def _alleles_from_contig(contig):
    '''Given a contig it returns a dict with the alleles as keys and the
    number of times they appear as values.
    It filters the alleles with less that min_num_of_reads.
    '''
    alleles = {}
    for row_index, sequence in enumerate(contig):
        if  sequence is None:
            continue 
        try:
            # This is to used when he sequence is a seqRecord
            # or a SeqWithQuality
            allele = sequence.seq.upper()
        except AttributeError:
            #This is to use when the sequence is a locatable sequence 
            allele = sequence.sequence.seq.upper()
            
        if allele not in alleles:
            alleles[allele] = []
        alleles[allele].append(row_index)
    #we filter the alleles that have not been read enough times
    for allele in alleles.keys():
        if len(alleles[allele]) < CONFIG.min_num_of_reads:
            del alleles[allele]
    return alleles

def _alleles_dict_to_set(alleles):
    ''' With this function we return a set of all reads and another set of
     reads that finish with indel ''' 
    # We get a list of alleles with indel as the last character
    indel_reads = set()
    # All the reads readed at least CONFIG.min_num_reads times 
    all_reads  = set()
    for allele in alleles:
        all_reads = all_reads.union(alleles[allele])
        if allele[-1] == CONFIG.indel_char:
            indel_reads = indel_reads.union(alleles[allele])
    return indel_reads, all_reads

def _build_location_from_index(index):
    ''' It returns alocation or a int depending on the length between
    start and end.
    The index can be an int or a Location.
    '''
    if isinstance(index, int):
        start = index
        stop = index + 1
    else:
        start = index.start
        stop = index.stop

    if start == stop - 1:
        loc = start
    else:
        loc = NonStaticParentLocation(start=start, end=stop - 1)
    return loc

def _seqvariation_in_subcontig(contig, colindex, rows):
    '''Given a contig a column index and a list of rows (seqs) it looks if
    there is a SeqVariation. If it founds one it returns it, otherwise it
    returns None.'''
    #which are the alleles in this subcontig?
    seqs      = [contig[seqindex, colindex] for seqindex in rows]
    subcontig = Contig(sequences=seqs)
    alleles   = _alleles_from_contig(subcontig)

    #do we have more than one allele == is a seqvar?
    if len(alleles) > 1:
        loc = _build_location_from_index(colindex)
        return SeqVariation(alleles=alleles, location=loc, alignment=contig)
    else:
        return None


def seqvariations_in_alignment(alignment):
    ''' We use this method to yield the SequenceVariation found in an alignment.
    
    The alignment (contig) must be a list of SeqRecord-like class objects.'''
    inchar     = CONFIG.indel_char
    alignment.return_empty_seq = True

    #we go through every column in the alignment
    col_index = 0
    while col_index < alignment.ncols:
        #which are the alleles in the column col_index?
        cseq = alignment[:, col_index: col_index + 1]
        # If all the secuences that we get in this column are none,
        # we have finished
        alleles = _alleles_from_contig(cseq)
        #are we dealing with an indel or with an snp?
        if inchar not in alleles:
            #the snp case is simple, we return it and we go for the next column
            if len(alleles) > 1:
                yield SeqVariation(alleles=alleles, location=col_index,
                                   alignment=alignment)
            col_index += 1
        else:
            # the indel case is complex because an indel can cover several
            # columns
            # We need to know where the indel finish, so we go to the following
            # columns looking for the end.
            indel_reads      = set(alleles[inchar])
            previous_alleles = alleles
            #we know look for the indel span
            col_indel_end = col_index + 1 #the last column covered by the indel
            while True:
                #if there are more columns
                if col_indel_end <= alignment.ncols:
                    #we need to know if the indel continues or not
                    #we calculate a set with the reads that are still streching
                    #the indel.
                    
                    #we need the alleles taking into account the whole span of
                    #the indel to deal with situations like:
                    #    --
                    #    -A
                    subcontig         = alignment[:, col_index:col_indel_end+1]
                    subcontig_alleles = _alleles_from_contig(subcontig)
                    #Which reads have an indel in the last column?
                    #subcontig_good_reads are the reads in this alignment
                    #        section that have an allele read at least 
                    #        min_num_reads times
                    #subcontig_indel_reads are the reads in the section that
                    #        we're considering now that finish with an indel
                    subcontig_indel_reads, subcontig_good_reads = \
                                         _alleles_dict_to_set(subcontig_alleles)
                    #So here we have the reads that are streching the indel
                    #These are the reads that finished with an indel in the
                    #previous iteration (the previous col_indel_end) and are
                    #still an indel
                    #        --T          --T
                    #        ---          TT-
                    #          ^- follows   ^- not follows
                    #pylint: disable-msg=C0301
                    cont_indels = indel_reads.intersection(subcontig_indel_reads)
                else:
                    # If there are no more columns. 
                    # So none of the indels  continues.  
                    cont_indels           = set()
                    # We need to save the last valid reads with indels.
                    subcontig_indel_reads = indel_reads
                    #We should take in account the alleles of the last iteration
                    subcontig_alleles     = previous_alleles
                if len(cont_indels) == 0: # Indel doesn't continue
                    # I do a new subcontig only with subcontig_indel_list aleles
                    # maybe we have to return a SeqVariation. It depends if 
                    # enough alleles have been read more than min_number_reads
                    # times 
                    colindex = slice(col_index, col_indel_end)
                    seqvar   = _seqvariation_in_subcontig(alignment, colindex, \
                                                 subcontig_good_reads)
                    if seqvar is not None:
                        yield seqvar
                        col_index = col_indel_end
                        # we have yield a seqvar, so we are done with the indel 
                        # streching
                        break
                    # To have to actualize the set of reads that still finish 
                    # with an indel in order to continue with the new iteration
                    indel_reads      = subcontig_indel_reads
                    previous_alleles = subcontig_alleles
                    # After elongating the indel we haven't found a seqvar
                    # So we are done with the indel streaching 
                    # Now We go to search for the sext seqvar in the alignment's
                    # next column
                    col_index += 1
                    break
                else:
                    # The indel is still growing
                    # We actualizes the list of reads that continues streching
                    # the indel
                    indel_reads    = subcontig_indel_reads
                    # We strech it one more column    
                    col_indel_end += 1
