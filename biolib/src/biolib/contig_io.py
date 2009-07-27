'Module to write out and to read in Contigs (alignments)'

from re import match
from biolib.contig import Contig,  locate_sequence
from biolib.seqs import SeqWithQuality, Seq
from biolib.biolib_utils import fasta_str, FileIndex
from Bio.Seq import Seq

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

def get_parser(fhand, format):
    'Given a file and a format it returns a suitable Contig parser'
    available_parsers = {'caf': CafParser, 'ace': AceParser,
                         'bowtie': BowtieParser}
    parser = available_parsers[format](fhand)
    return parser

def get_parser_by_name(fname):
    ''''Given a file returns suitable contig parser. Looking in file name'''
    if fname[-3:].lower() == 'ace':
        return AceParser(open(fname, 'r'))
    elif fname[-3:].lower() == 'caf':
        return  CafParser(open(fname, 'r'))

def contig_to_fasta(infpath, contig_list=None):
    '''It returns a fasta file content with the consensus sequences of a
    caf or ace file   '''
    fasta = ''
    parser   = get_parser_by_name(infpath)
    for contig in parser.contigs():
        name = contig.consensus.sequence.name
        if contig_list is  None or name in contig_list:
            sequence = contig.consensus.sequence
            fasta   += fasta_str(sequence, name)
    return fasta

def _pairwise_alignments(contig):
    '''It yields all the pairwise alignments present in the contig

    A contig is a multiple alignment.
    If there are only two sequences it will return the contig itself.
    '''
    if len(contig) == 2:
        yield contig

def _build_contig_from_dict(reads):
    'Given a dict with the contig info it returns a Contig'
    #now we can create the contig
    contig = Contig()
    for read_name in reads:
        read = reads[read_name]
        seq = Seq(read['dna'])
        start = read['start']
        real_name = read['name'] #different for the consensus
        if read_name == 'consensus':
            #consensus
            seq = SeqWithQuality(name=real_name, seq=seq,
                                 qual=read['quality'])
            con = \
              locate_sequence(seq, location=start, parent=contig)
            contig.consensus = con
        else:
            #reads
            if 'quality' in read:
                qual = read['quality']
            else:
                qual = None
            seq = SeqWithQuality(name=real_name, seq=seq, qual=qual)
            mask = read['mask']
            forward = read['forward']
            if forward:
                strand = 1
            else:
                strand = -1
            contig.append_to_location(seq, start=start, mask=mask,
                                      strand=strand, forward=forward)
    return contig

class CafParser(object):
    ''' This class is used to parse caf files.'''
    def __init__(self, fhand):

        '''The initialitation.
        keyword arguments:
            fhand : caf file
        '''
        self._fhand     = fhand
        self._qual_index = {}
        self._seq_index = {}
        self._dna_index = {}
        self._type_index = {}
        self._build_index()

    def contig_names(self):
        'It yields the contig names.'
        for name, kind in self._type_index.items():
            if kind == 'Is_contig':
                yield name

    def read_names(self):
        'It yields the read names.'
        for name, kind in self._type_index.items():
            if kind == 'Is_read':
                yield name

    def _build_index(self):
        '''It takes a caf file and after reading it it returns an index with
           section positions. We define section as the paragraph of text where
            each caf data type is represented . example:
                    DNA: name
                    aaaaaaaaaaaaaaa
                    ttttttttttt
                    cccccccccc
         It stores as well if the secuence is a contig or a read
         '''
        rawline  = "Filled"
        sec_in = False
        fhandler = self._fhand
        while len(rawline) != 0:

            prior_tell = fhandler.tell()
            rawline    = fhandler.readline()
            line       = rawline.strip()

            #if match("\w*\s*:\s*\w*", line):
            mode = rawline.split(':', 1)[0].strip().lower()
            if mode in ('dna', 'basequality', 'sequence'):
                #pylint: disable-msg=W0612
                mode_, name = line.split(":")
                name = name.strip()
                if mode == "sequence":
                    sec_in = True
                else:
                    sec_in = False

                if mode == 'dna':
                    self._dna_index[name] = prior_tell
                elif mode == 'basequality':
                    self._qual_index[name] = prior_tell
                elif mode == 'sequence':
                    self._seq_index[name] = prior_tell

            if sec_in:
                if line == "Is_read":
                    self._type_index[name] = "Is_read"
                elif line == "Is_contig":
                    self._type_index[name] = "Is_contig"

    def contigs(self):
        '''It returns a generator that yields the contigs.'''

        for name in self._seq_index:
            if self._type_index[name] == 'Is_contig':
                yield self.contig(name)

    def reads(self):
        '''It returns a generator with the reads'''

        for seq_rec_name in self._seq_index:
            if self._type_index[seq_rec_name] == 'Is_read':
                yield  self.read(seq_rec_name)

    def _return_section(self, position):
        ''' It returns a section giving a position in the file. It will take
        the text until it finds the next statement ( DNA: , Sequence: or
        BaseQuality '''
        content = []
        fhandler = self._fhand
        fhandler.seek(position)

        line = fhandler.readline()# To pass the header
        line = fhandler.readline()
        while True:
            line = line.strip()
#            if (len(line) == 0) or (match("\w*\s*:\s*\w*", line)):
            if (len(line) == 0) or \
            (match("^[DNA|Sequence|BaseQuality|BasePosition]\s*:", line)):
                break
            content.append(line)
            line = fhandler.readline()
        return content

    def _read_dna_section(self, sec_rec_name):
        ''' It returns the dna secuence in a string.It needs the sec_rec name
        '''
        dna_pos     = self._dna_index[sec_rec_name]
        dna_section = self._return_section(dna_pos)
        dna = ''
        for line in dna_section:
            dna += line.strip()
        return dna

    def _read_qual_section(self, sec_rec_name):
        ''' It returns the base quality list. It needs the sec_rec name.'''
        try:
            base_quality_pos = self._qual_index[sec_rec_name]
        except KeyError:
            raise ValueError('No quality for the given read')
        base_quality_section = self._return_section(base_quality_pos)

        base_quality = []
        for line in base_quality_section:
            line          = line.strip()
            base_quality += line.split(' ')
        return base_quality

    @staticmethod
    def _get_align_to_scf(line):
        ''' It reads Alig_to_SCF line and returns a tupla with the info
         structured. Tupla order: scf_start, scf_end, read_start, read_end '''
        item       = line.split(" ")
        return (int(item[1]), int(item[2]), int(item[3]), int(item[4]))

    @staticmethod
    def _get_assembled_from(line):
        ''' It reads Assembled_from line and returns 2 elements: The first is
        the read name, and the second one is a dict with 4 keys:
        c1, c2, r1, r2'''
        item         = line.split(" ")
        #pylint: disable-msg=C0103
        #Assembled_from Read_name c1 c2 r1 r2
        name    = item[1]
        c1      = int(item[2])
        c2      = int(item[3])
        r1      = int(item[4])
        r2      = int(item[5])
        coord_dict = {'c1':c1, 'c2':c2, 'r1':r1, 'r2':r2}
        return name, coord_dict

    @staticmethod
    def _calc_read_coords_in_contig(coord_dict, read_length, read_info):
        '''Given the coords from the assembled_from line (c1, c2, r1, r2) it
        calculates where the read starts in the contig, if the read is forward
        or not and if the masked region in the read.
        read_info is a dict and the results will be added there.
        '''
        #we need to know where the read starts in the contig and its mask
        #pylint: disable-msg=C0103
        r1 = coord_dict['r1']
        r2 = coord_dict['r2']
        c1 = coord_dict['c1']
        c2 = coord_dict['c2']
        mask = (r1 - 1, r2 - 1)
        forward = True
        if c1 < c2:
            contig_start = c1 - r1 + 1
        else:
            forward = False
            contig_start = c2 - (read_length - r2)
        #contig start
        read_info['start'] = contig_start
        read_info['mask'] = mask
        read_info['forward'] = forward

    def _read_seq_section(self, sec_rec_name):
        ''' It return a dictionary with the content of the secuence section.

        It stores the key and the value, but there are  2 exceptions:
        1) reads key's value is a dict with the reads. Where the key is
        the read name and the value are a list of tuples where each tupla
        contains:
            (contig_start, contig_end, read_start,read_end)

        2) scf_alignments key's value is a dict with the alignemnts.
         Where the key is the SCF file name and the value is a list of
         tuples where each tupla contains:
            (scf_start, scf_end, read_start, read_end)
        '''

        # variable type in this section. All of them are in array to easy use
        # of them
        seq_type = ('Is_read', 'Is_contig', 'Is_group', 'Is_assembly')
        state    = ('Padded', 'Unpadded')
        sec_info = {}
        reads    = {}
        scf_alignment = []

        sequence_pos     = self._seq_index[sec_rec_name]
        sequence_section = self._return_section(sequence_pos)

        for line in sequence_section:
            line = line.strip()
            if line in state:
                sec_info['state'] = line
            elif line in seq_type:
                sec_info['type'] = line
            elif line.startswith('Align_to_SCF'):
                scf_alignment.append(self._get_align_to_scf(line))
            elif line.startswith('Assembled_from'):
                read_name, coords = self._get_assembled_from(line)
                if read_name not in reads:
                    reads[read_name] = []
                reads[read_name].append(coords)
            else:
                items = line.split(' ')
                sec_info[items[0]] = " ".join(items[1:])

        sec_info['name'] = sec_rec_name
        if reads:
            sec_info['reads'] = reads
        else:
            scf_file_name     = sec_info['SCF_File']
            scf_alignments = {}
            scf_alignments[scf_file_name] = scf_alignment
            sec_info['scf_alignments'] = scf_alignments
        return sec_info

    def read(self, name):
        '''Given a read name it returns a SeqWithQuality object'''

        if self._type_index[name] == 'Is_read':

            seq_info = self._read_seq_section(name)
            dna      = self._read_dna_section(name)
            quality  = self._read_qual_section(name)

            seq_rec  = SeqWithQuality(seq=Seq(dna), name=name, qual=quality )
            seq_rec.annotations =  seq_info
            return seq_rec
        else:
            raise RuntimeError (name + 'does not corresponds to a read')

    @staticmethod
    def _contig_start(reads):
        ''' It calculates which is the position in which the most left read
        start.
        '''
        contig_start = None
        for read in reads:
            contig_start_r = reads[read]['start']
            if contig_start is None:
                contig_start = contig_start_r
            elif contig_start > contig_start_r:
                contig_start = contig_start_r
        return contig_start

    @staticmethod
    def _move_reads(reads, amount_to_move):
        'It moves the reads the given a amount of bp to the right.'
        for read in reads:
            reads[read]['start'] += amount_to_move

    @staticmethod
    def _read_mask(read_lenght, annot_dict):
        ''' It returns the mask of the read. It gets the information using
         the annotations dictionary'''
        read = range(1, read_lenght)
        for key in annot_dict:
            if key.lower() == 'clippping' or key.lower() == 'seq_vec':
                start = annot_dict[key].split(" ")[1]
                end   = annot_dict[key].split(" ")[2]
                for i in range(int(start), int(end)):
                    if i in read:
                        read.remove(i)
        mask_start = read[0]
        mask_end   = read[-1]
        return [mask_start, mask_end]

    def contig(self, name):
        '''Given a name it returns a Contig'''
        #is the given name a contig?
        try:
            kind = self._type_index[name]
        except KeyError:
            raise ValueError('Given name is not a contig')
        if kind != 'Is_contig':
            raise ValueError('Given name is not a contig')

        #to store the seq and qual from the reads and the consensus
        reads = {}
        #we read the file and we get a dict for the consensus
        contig_info               = self._read_seq_section(name)
        reads['consensus']        = {}
        reads['consensus']['dna'] = self._read_dna_section(name)
        try:
            qual = self._read_qual_section(name)
        except ValueError:
            qual = None
        reads['consensus']['quality']    = qual
        reads['consensus']['start']   = 1 #the consense always starts at 1
        reads['consensus']['mask']    = None
        reads['consensus']['forward'] = True
        reads['consensus']['name']    = name
        for read_name in contig_info['reads']:
            #we need the dna seq and the quality
            dna = self._read_dna_section(read_name)
            try:
                qual = self._read_qual_section(read_name)
            except ValueError:
                qual = None
            #we store all that info in a dict
            reads[read_name]         = {}
            reads[read_name]['name'] = read_name
            reads[read_name]['dna']  = dna
            reads[read_name]['quality'] = qual
            #where does every sequence start in the contig and
            #which is it's mask?
            assembled_from_line = contig_info['reads'][read_name]
            if len(assembled_from_line) != 1:
                msg = 'We do not support unpadded sequences'
                raise NotImplementedError(msg)
            read_len = len(dna)
            self._calc_read_coords_in_contig(assembled_from_line[0], read_len,
                                             reads[read_name])
        #some reads might be negative, we want the most negative one moved
        #to the position 0 because our Contig does not support negative values
        contig_start = self._contig_start(reads)
        #now we move the read, because the most left one should be at 0
        self._move_reads(reads, -contig_start)
        return _build_contig_from_dict(reads)

class AceParser(object):
    '''It parses the ace assembly files.'''
    def __init__(self, fhand):
        '''It initializes the ace parser with an ace file name.'''
        self._fileh = fhand
        self._contig_index = {} #in which byte is every contig located
        self._contig_names = []      #a list with the contig names
        self._read_names   = []      #a list with the read names
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
                self._contig_names.append(contig_name)
            if line.startswith('RD') and line[2].isspace():
                read_name = line.split()[1]
                self._read_names.append(read_name)
        self._contig_index = index

    def contig_names(self):
        'It yields the contig names.'
        for name in self._contig_names:
            yield name

    def read_names(self):
        'It yields the read names.'
        for name in self._read_names:
            yield name

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
            if base == '-':
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
        dna = ''.join(co_items[6:])
        #we have to change the * for -
        dna = dna.replace('*', '-')
        reads['consensus']['dna'] = dna
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
            #we have to replace the * with -
            dna = ''.join(items[5:])
            dna = dna.replace('*', '-')
            reads[name]['dna'] = dna
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
        return reads


    def contig(self, name):
        '''Given a contig name it returns a Contig instance'''
        #we go to the location where the Contig is in the file
        try:
            file_position = self._contig_index[name]
        except KeyError:
            raise ValueError('Contig not in file: ' + name)
        #we're in this contig until we find the next one or the end of the file
        contig_in_file = self._read_contig_from_file(file_position)
        # Now we convert the info parsed to our standar dict of saving contigs
        contig_dict    = self._build_contig(contig_in_file)
        return _build_contig_from_dict(contig_dict)

    def contigs(self):
        '''It yields all contigs.'''
        for name in self._contig_index:
            yield self.contig(name)

class BowtieParser(object):
    'It parsers a bowtie map file and it prepares the contigs from it'
    def __init__(self, fhand):
        'It requires a bowtie.map file to create the contigs'
        self._fhand = fhand
        self._index = FileIndex(fhand,
                                item_start_patterns=['^'],
                                key_patterns=['^([^\t]+)\t'],
                                type_patterns=['\t[+|-]\t([^ \t]+)'])

    @staticmethod
    def _parse_read(read_line):
        '''Given a bowtie line it returns a read and its contig location'''
        #pylint: disable-msg=W0612
        (read_name, orientation, contig_name, read_location, sequence,
                    qualities) = read_line.split('\t')[:6]
        read_location = int(read_location)
        # The encoded quality values are on the Phred scale and the encoding is
        #ASCII-offset by 33 (ASCII char !
        qualities = [ord(qual) - 33 for qual in qualities]
        read = SeqWithQuality(name=read_name, seq=sequence, qual=qualities)
        return read, read_location

    def contig(self, name):
        'Given a contig name it returns the contig'
        contig = Contig()
        for read in self._index[name].items():
            read_content = self._index[name][read]
            read, location = self._parse_read(read_content)
            contig.append_to_location(sequence=read, start=location)
        return contig

    def contigs(self):
        'It yields all contigs'
        for name in self._index.types():
            yield self.contig(name)
