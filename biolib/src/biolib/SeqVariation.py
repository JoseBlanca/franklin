'''
Created on 2009 mar 25

@author: peio
'''
#from biolib.contig import Contig, Location

class SeqVarConf(object):
    '''
    This class contains some switches to configure to your needs
    '''
    def __init__(self, min_num_of_reads = 1, only_snp = False):
        ''' Here we initialize teh object with the configuration we want'''
        self._min_num_of_reads = min_num_of_reads
        self._only_snp         = only_snp
        
    def min_num_reads(self):
        ''' Returns the minimun number of reads for the second most
        frecuent alele'''
        return self._min_num_of_reads
    
    def only_snp(self):
        ''' Returns if we are going to use only snps or not'''
        return self._only_snp 
        
class SeqVariation(object):
    '''
    This class is used to represent any kind of sequence variation in 
    and alignment. The variation can be a Snp or a InDel, or ...
    '''

    def __init__(self, name, aleles, location = None, alignment = None):
        '''
        This class collect all the information about each secuence variation
        in an alignment.Sequence variations could be SNPs or InDels. There is a
        filter in the conf class to use only SNPs in this class.
        Arguments:
            name     - Name of the variations
            location - Location of the Seq variation. Location class object
            alignemt - Alignment where it procedes. It can be a COntig class 
                       object
            Aleles - . The aleles contains a dict with the reads and the aleles
                       it have. example:
                        aleles= {read1 = 'A', read2 = 'T',...}
                        Keys are read names
                        Value could be a string or a biopython Seq class object 
        To initialize you only need the aleles information and the name.
        '''
        
        self._name      = name
        self._aleles    = aleles
        self._location  = location
        self._alignment = alignment
    
    
    
    def __repr__(self):
        ''' It prints a readable SeqVariation information '''
        
        pstring  = 'SequenceVariation: ' + str(self._name) + '\n'
        if self._location:
            pstring += self._location.__repr__() + '\n'
        if self._alignment:
            pstring += self._alignment.__repr__() + '\n'
        pstring += self._aleles.__repr__()
        return pstring
    
    def aleles(self):
        ''' It returns a dict with the aleles and  the reads it contains'''
        return self._aleles
    
    def alignment(self):
        '''It returns the contig it belong to, The reference '''
        return self._alignment
        
    def name(self):
        '''It returns the name of the SeqVariation '''
        return self._name
    
def seq_var_in_alignment(contig):
    ''' We use this method to catch the Sequence variation from an
     alignment. The alignment (contig) MUST BE a list of Biopython SeqRecord
     class objects'''      
    
    ncols     = contig.ncols
    
    for loc_order in range(ncols):
        colum_seq = contig[:,loc_order]
        yield  colum_seq
         
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    