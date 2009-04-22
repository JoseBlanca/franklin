'''It parses the ace assembly files.'''

from biolib.contig import Contig, locate_sequence
from biolib.Seqs import SeqWithQuality
from Bio.Seq import Seq

class AceParser(object):
    '''It parses the ace assembly files.'''
    def __init__(self, fname):
        '''It initializes the ace parser with an ace file name.'''
        self._fileh = open(fname)
        self._contig_index = {} #in which byte is every contig located
        self._build_index()

    def _build_index(self):
        '''It scans the ace file and it creates an index.
        
        In the index the file position in which every contig starts is located.
        '''
        fileh = self._fileh
        fileh.seek(0)
        index = {}
        line = ' '
        while len(line) != 0:
            file_position = fileh.tell()
            line = fileh.readline()
            if line.startswith('CO') and line[2].isspace():
                contig_name = line.split()[1]
                index[contig_name] = file_position
        self._contig_index = index
 
    def _read_contig_from_file(self, position):
        '''Given a position in the file it reads the next contig.
        
        The position should be the first line of a contig in the file
        It returns a dict with the sections as keys and the section contents
        as values.
        '''
        #CO 1 30502 510 273 U
        #CCTCTCC*GTAGAGTTCAACCGAAGCCGGTAGAGTTTTATCACCCCTCCC
        #
        #BQ
        #20 20 20 20 20 20 20 20 20 20 20 20 20
        #
        #AF TBEOG48.y1 C 1
        #
        #BS 1 137 TBEOG48.y1
        #
        #RD TBEOG48.y1 619 0 0
        #CCTCTCC*GTAGAGTTCAACCGAAGCCGGTAGAGTTTTATCACCCCTCCC
        #
        #QA 1 619 1 619
        self._fileh.seek(position)
        line = ' '
        section = None  #file sections like contig (CO), read(RD), etc
        read_cache = {}  #it holds everything that has been read in this
                         #section, keys are the sections
        read = None     #name of the sequence read we're dealing
        while len(line) != 0:
            line = self._fileh.readline()
            #are we in a new section?
            if line[:2] in ('CO', 'BQ', 'AF', 'RD', 'QA') and line[2].isspace():
                section = line[:2]
                #Which read are we dealing with?
                if section in ('AF', 'RD'):
                    read = line.split()[1]
                elif section == 'BS':
                    read = line.split()[3]
                elif section in ('CO', 'BQ'):
                    read = 'consensus'
                #some of this sections can be found several times (one per read)
                #so they will be dicts
                if section not in read_cache:
                    read_cache[section] = {}
                read_cache[section][read] = ''
            #is this section already finished
            if line.isspace() or len(line) == 0:
                section = None
            #is this already another contig?
            if section != 'CO' and line.startswith('CO') and line[2].isspace():
                line = ''
            if section:
                read_cache[section][read] += line
        return read_cache

    @staticmethod
    def _fix_consensus_quality(sequence, quality):
        '''It adds a 0 to the quality list in the places where an * is found in
        the sequence.
        
        If at the end quality an seq have a different length it will raise a
        RuntimeError'''
        quality = quality[:] #we copy to leave the original intact
        for index, base in enumerate(sequence):
            if base == '*':
                quality.insert(index, 0)
        if len(sequence) != len(quality):
            msg = 'Consensus seq and qual do not match when qual is fixed.'
            raise RuntimeError(msg)
        return sequence, quality

    @staticmethod
    def _move_contig_dict(contig, num_bases):
        'It moves the contig dict to the right the given number of bases'
        for read in contig:
            contig[read]['start'] = contig[read]['start'] + num_bases
            #No read is suposed to start at a negative index
            if contig[read]['start'] < 0:
                msg = 'A read starts at a negative position after fix it'
                raise RuntimeError(msg)

    def _build_contig(self, contig_dict):
        '''Given a dict with the contig contents from an ace file it returns a
        Contig.'''
        #the format for CO is
        #CO name 30502 510 273 U \nsequence\nsequence
        #CO name, numbases, numreads, num_base_segments, complemented, sequence
        #the complemented letter means Uncomplemented, complemented
        co_items = contig_dict['CO']['consensus'].strip().split()
        reads = {}
        reads['consensus'] = {}
        reads['consensus']['name'] = co_items[1]
        reads['consensus']['complemented'] = co_items[5]
        reads['consensus']['dna'] = ''.join(co_items[6:])
        #In the ace format the consensus always start at the position 1, (0 in
        #the 0 based python indexing scheme)
        reads['consensus']['start'] = 0
        #BQ encodes the consensus quality
        reads['consensus']['quality'] = \
                                    contig_dict['BQ']['consensus'].split()[1:]
        #fix padded qualities
        #Consensus quality values are provided for the bases alone, the gaps
        #not being represented
        seq = reads['consensus']['dna']
        qual = reads['consensus']['quality']
        seq, qual = self._fix_consensus_quality(seq, qual)
        reads['consensus']['quality'] = qual

        #which are the reads?
        read_names = contig_dict['AF'].keys()
        #contig_start tracks where the contig start (index most to the left)
        contig_start = None
        for name in read_names:
            #AF TBEOG48.y1 C 1
            #The AF lines (one per aligned read) contain information of whether
            #the read is complemented (C) or not (U) followed by a 1-based
            #offset in the consensus sequence.  Note that the offset refers to
            #the beginning of the entire read in the alignment, not just the
            #clear range.  Thus the read acaggATTGA will have an offset of 1
            #even though the consensus truly starts at position 6.
            items = contig_dict['AF'][name].split()
            reads[name] = {}
            reads[name]['name'] = name  #just to be compatible with the
                                        #consensus dict structure
            reads[name]['complemented'] = items[2]
            #in fact all reads are forward now
            reads[name]['forward'] = True
            #we count from 0 so: start = start - 1
            reads[name]['start'] = int(items[3]) - 1
            #which is the minimum start
            if contig_start is None:
                contig_start = reads[name]['start']
            if contig_start > reads[name]['start']:
                contig_start = reads[name]['start']
            #the read sequence
            items = contig_dict['RD'][name].split()
            reads[name]['dna'] = ''.join(items[5:])
            #QA, quality clip
            #QA line contains two 1-based ranges. The second range represents
            #the clear range of the read, with respect to the read sequence
            #(padded and potentially complemented) as provided in
            #the RD record
            items = contig_dict['QA'][name].split()
            #we count from zero, so the seqs indexes are minus 1
            reads[name]['mask'] = int(items[3]) - 1, int(items[4]) - 1

        #at this point the contig might start at a negative location, that's 
        #not supported by our Contig, so we move everything to the right
        self._move_contig_dict(reads, -contig_start)
            
        #now we can create the contig
        contig = Contig()
        for read_name in reads:
            seq = Seq(reads[read_name]['dna'])
            start = reads[read_name]['start']
            name = reads[read_name]['name'] #different for the consensus
            if read_name == 'consensus':
                seq = SeqWithQuality(name=name, seq=seq,
                                     qual=reads[read_name]['quality'])
                con = \
                  locate_sequence(seq, location=start, parent=contig)
                contig.consensus = con
            else:
                seq = SeqWithQuality(name=name, seq=seq)
                mask = reads[read_name]['mask']
                contig.append_to_location(seq, start=start, mask=mask)
        return contig

    def contig(self, name):
        '''Given a contig name it returns a Contig instance'''
        #we go to the location where the Contig is in the file
        try:
            file_position = self._contig_index[name]
        except KeyError:
            raise ValueError('Contig not in file: ' + name)
        #we're in this contig until we find the next one or the end of the file
        contig_in_file = self._read_contig_from_file(file_position)
        return self._build_contig(contig_in_file)
        
    def contigs(self):
        '''It yields all contigs.'''
        for name in self._contig_index:
            yield self.contig(name)

