'''
Created on 2009 mar 25

@author: peio
'''
from biolib.contig import NonStaticParentLocation,Contig

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

def seqvariation_alleles_with_list(alleles, name=None, location=None,
                                   alignment=None):
    '''This is an alternative init for the SeqVariation.
    
    Here each value in the alleles dict is a list instead of an int with the
    count.
    '''
    int_alleles = {}
    for allele in alleles:
        int_alleles[allele] = len(alleles[allele])
    return SeqVariation(int_alleles, name, location, alignment)

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
            pstring += self.location.__repr__() + '\n'
        if self.alignment:
            pstring += self.alignment.__repr__() + '\n'
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
            if alleles[allele] < CONFIG.min_num_of_reads:
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
        allele = sequence.seq.upper()
        if allele not in alleles:
            alleles[allele] = []
        alleles[allele].append(row_index)
    #we filter the alleles that have not been read enough times
    for allele in alleles.keys():
        if len(alleles[allele]) < CONFIG.min_num_of_reads:
            del alleles[allele]
    return alleles
  
def seqvariations_in_alignment(contig):
    ''' We use this method to catch the Sequence variation from an
     alignment. The alignment (contig) must be a list of SeqRecord-like
     class objects'''      
    inchar     = CONFIG.indel_char
    col_number = 0
    while col_number < contig.ncols:
        cseq = contig[:, col_number: col_number + 1]
        #which are the alleles?
        alleles = _alleles_from_contig(cseq)
        if inchar not in alleles:
            if len(alleles) > 1:
                yield seqvariation_alleles_with_list(alleles=alleles,
                                                     location=col_number,
                                                     alignment=contig)
            col_number += 1
        else:
            # We are finding allele's length, And we follow continuous indels
            indel_reads      = set(alleles[inchar])
            previous_alleles = alleles
            right_col_number = col_number + 1
            while True:
                #we need to know if the indel continues in the next column or
                #not if there is no next column the indel does not continue
                if right_col_number <= contig.ncols:
                    #here we create a new subcontig with one more column
                    subcontig         = contig[:, col_number:right_col_number+1]
                    subcontig_alleles = _alleles_from_contig(subcontig)
                    # We get a list of alleles with indel as the last character
                    subcontig_indel_reads = set()
                    # All the reads readed at least CONFIG.min_num_reads times 
                    subcontig_good_reads  = set()
                    for allele in subcontig_alleles:
                        subcontig_good_reads = \
                            subcontig_good_reads.union(subcontig_alleles[allele])
                        if allele[-1] == inchar:
                            subcontig_indel_reads = \
                            subcontig_indel_reads.union(subcontig_alleles[allele])
                    # These are the read number with the indels we are elongating
                    cont_indels = indel_reads.intersection(subcontig_indel_reads)
                else:
                    cont_indels = set()
                    subcontig_indel_reads = indel_reads
                    subcontig_alleles = previous_alleles
                if len(cont_indels) == 0:
                    #Indels doesn't continue, so
                    #I do a new subcontig only with subcontig_indel_list aleles
                    colindex          = slice(col_number, right_col_number)
                    good_reads_contig = [contig[seqindex, colindex] \
                                        for seqindex in subcontig_good_reads]
                    final_contig      = Contig(sequences=good_reads_contig)
                    #we have to recalculate the previous alleles using the
                    #good_reads_contig
                    previous_alleles = _alleles_from_contig(final_contig)

                    #do we have more than one allele == is a seqvar?
                    if len(previous_alleles) > 1:
                        if col_number == right_col_number - 1:
                            loc = col_number
                        else:
                            loc = NonStaticParentLocation(start=col_number,
                                                      end=right_col_number -1)
                        yield SeqVariation(alleles  = previous_alleles, 
                                           location = loc, 
                                           alignment= contig)
                        col_number = right_col_number
                        #we have yield a seqvar, so we go to the main while
                        break
                    indel_reads      = subcontig_indel_reads
                    previous_alleles = subcontig_alleles
                    #after elongating the indel we haven't found a seqvar
                    #we don't return a seqvar so colnumber goes to the next
                    #column
                    col_number += 1
                    break
                else:
                    #the indel is still growing
                    indel_reads        = subcontig_indel_reads
                    right_col_number += 1
