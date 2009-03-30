'''
Created on 2009 mar 25

@author: peio
'''
#from biolib.contig import Location

class _SeqVarConf(object):
    '''This class contains some switches to configure to your needs
    '''
    #No public methods at all, only some properties.
    #we could use a namedtuple but it would be less clear.
    #pylint: disable-msg=R0903
    def __init__(self, min_num_of_reads=2, only_snp=False,
                indel_char='-'):
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
        self.indel_char      = indel_char


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
    
def seq_var_in_alignment(contig):
    ''' We use this method to catch the Sequence variation from an
     alignment. The alignment (contig) MUST BE a list of Biopython SeqRecord
     class objects'''      
    for loc_order in range(contig.ncols ):
        colum_seq = contig[:, loc_order:loc_order+1]
        print len( colum_seq)
         
