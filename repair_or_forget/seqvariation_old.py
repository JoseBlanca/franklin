'''
Created on 2009 mar 25

@author: peio
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

from biolib.contig import Contig
from biolib.locatable_sequence import NonStaticParentLocation
from biolib.biolib_seqio_utils import temp_fasta_file
from biolib.biolib_utils import get_start_end
from biolib.biolib_cmd_utils import call
from biolib.seqs import SeqRecord

class _SeqVarConf(object):
    '''This class contains some switches to configure to your needs
    '''
    #No public methods at all, only some properties.
    #we could use a namedtuple but it would be less clear.
    #pylint: disable-msg=R0903
    def __init__(self, min_num_of_reads=2, only_snp=False,
                indel_char='-', empty_char='', valid_alleles=None):
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
        #too many arguments
        #pylint: disable-msg=R0913
        self.min_num_of_reads = min_num_of_reads
        self.only_snp         = only_snp
        self.indel_char       = indel_char
        self.empty_char       = empty_char
        if valid_alleles is None:
            valid_alleles  = ['A', 'a', 'T', 't', 'C', 'c', 'G', 'g', \
                              self.indel_char]
        self.valid_alleles = valid_alleles
        #Encimes to use with remap
        self.common_enzymes = ['ecori', 'smai', 'bamhi', 'alui', 'bglii',
                               'sali', 'bgli', 'clai', 'bsteii', 'taqi',
                               'psti', 'pvuii', 'hindiii', 'ecorv', 'xbai',
                               'haeiii', 'xhoi', 'kpni', 'scai', 'banii',
                               'hinfi', 'drai', 'apai', 'asp718']
CONFIG = _SeqVarConf()

def allele_count(allele):
    'It returns how many times an allele has been read.'
    if isinstance(allele, int):
        allele_count_ = allele
    else:
        allele_count_ = len(allele)
    return allele_count_

SNP = 0
INSERTION = 1
DELETION = 2

class SeqVariation2(object):
    '''
    This class is used to represent any kind of sequence variation in
    and alignment. The variation can be a Snp or a InDel, or ...
    '''
    def __init__(self, alleles, reference, kinds, qualities=None, name=None,
                 location=None, ):
        self.alleles = alleles
        self.reference = reference
        if qualities is None:
            self.qualities = []
        else:
            self.qualities = qualities
        self.name = name
        self.location = location
        self.type = type
    def __str__(self):
        return '%s %d: %s:' % (self.reference, self.location, str(self.alleles))


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
            alleles - . The alleles contains a dict with the alleles and the
                        number of times that it has been read.
                        example:
                            alleles= {'A':1, 'T':4}
            name     - Name of the variation (default None)
            location - Location of the SeqVariation. Location class object or
                       int (default None).
            alignment - Alignment where it procedes. It can be a Contig class
                       object (default None).
        '''
        self.name         = name
        self._num_reads   = None  # the alleles dict with the alleles as keys
        self._set_alleles(alleles)
        self.location     = location
        self.alignment    = alignment
        self.annotations  = {}

    def __repr__(self):
        ''' It prints a readable SeqVariation information '''

        pstring  = 'SequenceVariation: ' + str(self.name) + '\n'
        if self.location:
            pstring += 'Location:' + self.location.__repr__() + '\n'
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
            allele_count_ = allele_count(alleles[allele])
            if allele_count_ < CONFIG.min_num_of_reads:
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

    def sorted_alleles(self):
        '''It returns a list of sorted alleles according to the number of times
        that have been read'''
        #first we build the list
        alleles = []
        alleles_dict = self.alleles
        for name, times in alleles_dict.items():
            alleles.append((name, times))
        #now we sort it
        return sorted(alleles,
                      lambda x, y:allele_count(y[1]) - allele_count(x[1]))

def seqvar_summary(seqvar):
    '''It takes a seqvar and summarizes, returning a library file like format
     string'''
    contig_name = seqvar.alignment.consensus.name
    sorted_alleles = seqvar.sorted_alleles()
    second_alleles = allele_count(sorted_alleles[1][1])

    kind               = seqvar.kind()
    loc_start, loc_end = get_start_end(seqvar.location)
    id_snp             = contig_name + '_' +  str(loc_start)
    if 'pic' not in seqvar.annotations:
        calculate_pic(seqvar)
    if 'enzyme' not in seqvar.annotations:
        cap_enzime(seqvar)
    pic                = seqvar.annotations['pic']
    enzymes            = seqvar.annotations['enzymes']

    toprint  = "snp\n"
    toprint += "\tname:%s\n" % id_snp
    toprint += "\tcontig:%s\n" % contig_name
    toprint += "\tstart:%d\n" % loc_start
    toprint += "\tend:%d\n" % loc_end
    toprint += "\tkind:%s\n" % kind
    toprint += "\tannotations: pic:%f\nsecond_allele:%d\nenzymes:%s\n" % \
                                                (pic, second_alleles, enzymes)
    alleles_toprint = ''
    for allele, reads in seqvar.alleles.items():
        reads = [seqvar.alignment[read].name for read in reads]
        alleles_toprint += "%s:%s" % (allele, ','.join(reads))
    toprint += "\talleles:%s\n" % alleles_toprint

    return toprint

def _alleles_from_contig(contig):
    '''Given a contig it returns a dict with the alleles as keys and the
    number of times they appear as values.
    It filters the alleles with less that min_num_of_reads.
    It filters the alleles we don't want
    '''
    alleles = {}

    for row_index, sequence in enumerate(contig):
        # This is to be able to deal with seqrecords or with strings
        # used in list of reads
        if  sequence is None:
            continue
        # We do it to support list of strings as a contig
        try:
            if sequence.isspace():
                continue
            #pylint: disable-msg=W0704
        except AttributeError:
            pass
        #Here we check if the allele have, valid nucleotides, defined in
        #the configuration
        have_novalid_nucleotide = False
        for letter in sequence:
            if str(letter) not in CONFIG.valid_alleles:
                have_novalid_nucleotide = True
        if have_novalid_nucleotide:
            continue

        allele = str(sequence).upper()
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
    ''' It returns a location or a int depending on the length between
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

def _seqvariation_in_subcontig(contig, alignment, colindex, rows):
    '''Given a contig a column index and a list of rows (seqs) it looks if
    there is a SeqVariation. If it founds one it returns it, otherwise it
    returns None.'''
    #which are the alleles in this subcontig?
    seqs = []
    for seqindex in rows:
        seqs.append(contig[seqindex][colindex])

    subcontig = Contig(sequences=seqs)
    alleles   = _alleles_from_contig(subcontig)

    #do we have more than one allele == is a seqvar?
    if len(alleles) > 1:
        loc = _build_location_from_index(colindex)
        return SeqVariation(alleles=alleles, location=loc, alignment=alignment)
    else:
        return None

def _select_colum_from_list(alignment, col_start, col_end):
    ''' Given a list and , row start and row end it returns a list with
    items from each row in given column(s)'''
    new_alignment = []
    for read in alignment:
        item = read[col_start:col_end]
        if not item:
            item = ' ' * (col_end - col_start)
        new_alignment.append(item)
    return new_alignment

def longest_read(alignment):
    ''' It returns the longest string lenght in the list'''
    longest = 0
    for read in alignment:
        len_read = len(read)
        if len_read > longest:
            longest = len_read
    return longest

def seqvariations_in_alignment(alignment):
    ''' We use this method to yield the SequenceVariation found in an alignment.

    The alignment (contig) must be a list of SeqRecord-like class objects.'''
    #the proxycontig strategy is used because is much much faster than
    #using our contig object directly
    proxycontig = contig_to_read_list(alignment)
    inchar     = CONFIG.indel_char
    try:
        alignment.return_empty_seq = True
        #pylint: disable-msg=W0704
    except AttributeError:
        pass
    ncols = alignment.ncols

    #we go through every column in the alignment
    col_index = 0
    while col_index < ncols:
        #which are the alleles in the column col_index?
        cseq = _select_colum_from_list(proxycontig, col_index, col_index + 1)

        # If all the secuences that we get in this column are none,
        # we have finished
        alleles = _alleles_from_contig(cseq)
        #are we dealing with an indel or with an snp?
        if len(alleles) < 2:
            col_index += 1
            continue
        if inchar not in alleles:
            #the snp case is simple, we return it and we go for the next column
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
                if col_indel_end <= ncols:
                    #we need to know if the indel continues or not
                    #we calculate a set with the reads that are still streching
                    #the indel.

                    #we need the alleles taking into account the whole span of
                    #the indel to deal with situations like:
                    #    --
                    #    -A
                    subcontig = _select_colum_from_list(proxycontig,
                                                            col_index,
                                                            col_indel_end + 1)

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
                    seqvar = _seqvariation_in_subcontig(proxycontig,
                                                        alignment,
                                                        colindex,
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

def seqvars_in_contigs(contigs_iter):
    'It returns seqvars found in an contig iterator'
    for contig in contigs_iter:
        for seqvar in seqvariations_in_alignment(contig):
            yield seqvar

def contig_to_read_list(contig):
    ''' It takes a contig class object and it fill a list with the reads.
    All the reads are '''
    reads = []
    for read in contig:
        reads.append(str(read))
    return reads


#functions to characterize the sequence variations
def calculate_pic(seq_variation):
    '''It calculates and returns the Polymorphic Information Content.

    The PIC was defined in Botstein 1980. Am. J. Hum. Genet. 32, 314 331 as  the
    probability that a given marker genotype of an offspring of an affected
    parent will allow deduction of the parental genotype at the marker locus.
    The calculation is done following Shete et al. 2000 Theoretical Population
    Biology 57, 265 271.
    '''
    # If it is already calculated it returns it
    if 'pic' in seq_variation.annotations:
        return seq_variation.annotations['pic']

    #pylint: disable-msg=C0103
    def _pic_sum_1(alleles, num_reads, num_alleles):
        '''It returns the first summation for the pic calculation'''
        # P is the frecuency one allele have been read
        suma  = 0
        for i in range(num_alleles):
            frec = (alleles[i] / num_reads) ** 2
            suma += frec
        return suma

    def _pic_sum_2(alleles, num_reads, num_alleles):
        '''It returns the second summation for the pic calculation '''
        suma = 0
        for i in range(num_alleles - 1):
            freci = (alleles[i] / num_reads) ** 2
            for j in range(i + 1, num_alleles):
                frecj = (alleles[j] / num_reads) ** 2
                suma += freci * frecj
        return suma


    alleles = seq_variation.alleles
    #the alleles can have the count or a list with the alleles, we make sure
    #that all have a count, and we convert the dict to a list
    alleles = [float(allele_count(allele)) for allele in alleles.values()]
    #how many reads are in total?
    num_reads = float(sum(alleles))
    #how many alleles are in total
    num_alleles = len(alleles)
    sum_1 = _pic_sum_1(alleles, num_reads, num_alleles)
    sum_2 = _pic_sum_2(alleles, num_reads, num_alleles)


    pic = 1.0 - sum_1 - ( 2 * sum_2)
    # We add the pik to the snp annotations
    seq_variation.annotations['pic'] = pic
    return pic

def cap_enzime(snp, all_enzymes=False):
    ''' It looks in the 2 most frecuent alleles if there is each of the enzimes
    cut diferently'''
    # If it is already calculated we return them
    if 'enzymes' in snp.annotations:
        return snp.annotations['enzymes']


    location = snp.location
    loc_start, loc_end = get_start_end(location)
    # It takes the two most frecuent alleles
    alleles = snp.sorted_alleles()
    allele1_orig = alleles[0][0]
    allele2_orig = alleles[1][0]
    inchar = CONFIG.indel_char
    if inchar in allele1_orig:
        allele1 = allele1_orig.replace(inchar, '')
    else:
        allele1 = allele1_orig
    if inchar in allele2_orig:
        allele2 = allele2_orig.replace(inchar, '')
    else:
        allele2 = allele2_orig

    # Now we need to know with sequence piece take from the consensus
    # How big is the piece to use with remap?
    piece_from_location = 7
    piece_start = loc_start - piece_from_location
    if piece_start < 0 :
        return None
        #raise ValueError(' The snp is too close to begining of consensus')
    piece_end = loc_end + piece_from_location
    consensus = str(snp.alignment.consensus)
    #the base sequence
    seq1      = consensus[piece_start: piece_end + 1]
    #the allele1 and 2 sequences
    al_start = loc_start - piece_start
    al_stop  = al_start + len(allele1)
    seq1 = seq1[:al_start] + allele1 + seq1[al_stop:]
    seq2 = seq1[:al_start] + allele2 + seq1[al_stop:]
    if len(seq1.strip()) < (piece_from_location * 2 + 1):
        return None
        #raise ValueError('The snp is in the end of the consensus')
    seq1 = SeqRecord(name='seq1', seq=seq1)
    seq2 = SeqRecord(name='seq2', seq=seq2)

    enzymes1 = _remap_run(seq1, all_enzymes)
    enzymes2 = _remap_run(seq2, all_enzymes)

    enzymes = set(enzymes1).symmetric_difference(set(enzymes2))

    # We add the enzymes to the snp annotations
    snp.annotations['enzymes'] = enzymes
    return list(enzymes)


def _remap_run(seq, all_enzymes):
    '''this command runs remap EMBOSS binary and returns ...'''
    # Minimun length of the restriction enzyme recognition site
    sitelen = 4
    seq_file = temp_fasta_file(seq)
    seq_filename = seq_file.name

    if all_enzymes:
        enzymes = 'all'
    else:
        enzymes = ",".join(CONFIG.common_enzymes)

    cmd = ['remap', '-sequence', seq_filename, '-enzymes', enzymes,
           '-sitelen' , str(sitelen), 'stdout']

    try:
        stdout, stderr, retcode  = call(cmd)
    except OSError:
        raise OSError('remap binary does not exits or it is not in path')

    if retcode:
        raise RuntimeError('remap error: '+ stderr)

    return  _parse_remap_output(stdout)

def _parse_remap_output(remap_output):
    ''' It takes the remap output and it returns a set list with the enzymes
     that cut there'''
    section = ''
    enzymes = []
    for line in remap_output.split('\n'):
        line = line.strip()
        if line.isspace() or len(line) < 2:
            continue
        if section == 'cut':
            if line.startswith('#'):
                section = ''
            else:
                enzymes.append(line.split()[0])

        if line.startswith('# Enzymes that cut'):
            section = 'cut'
            continue

    return enzymes



