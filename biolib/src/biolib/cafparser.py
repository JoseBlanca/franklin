'''
Created on 2009 mar 9

@author: peio

Script who reads a caf file and takes the information for each
contig. It creates a file for each contig to  easyly use them
'''

from re import match
from biolib.contig import Contig, locate_sequence
from biolib.Seqs import SeqWithQuality
from Bio.Seq import Seq

class CafParser(object):
    ''' This class is used to parse caf files.'''
    def __init__(self, fname):

        '''The initialitation.
        keyword arguments:
            fname : caf file name
        '''
        self._fname     = fname
        self._qual_index = {}
        self._seq_index = {}
        self._dna_index = {}
        self._type_index = {}
        self._build_index()

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
        
        fhandler = open(self._fname,'rt')
        rawline  = "Filled"
        sec_in = False

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
        
        for seq_rec_name in self._seq_index:
            if self._type_index[seq_rec_name] == 'Is_contig':
                yield self.contig(seq_rec_name)
    
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
        fhandler = open(self._fname, 'r')
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
            contig_start = c2 - (read_length - r2) - 1
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

    def _get_seq_rec_full(self, seq_rec_name):
        ''' It returns the complete info of a sec_record. It uses the index 
         to access the file so we do no t need to read the whole file.
         We need the name of the sec record '''
        
        seq_rec_info = self._read_seq_section(seq_rec_name) 
        # First we take dna secuence
        seq_rec_info['DNA'] = self._read_dna_section(seq_rec_name)
        # If BaseQuality is in the sec record
        if seq_rec_name in self._qual_index:
            seq_rec_info['BaseQuality'] = self._read_qual_section(seq_rec_name)
        
        return seq_rec_info
    
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
        reads['consensus']['qual']    = qual
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
            reads[read_name]['qual'] = qual
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
        print contig_start
        self._move_reads(reads, -contig_start)
        return self._build_contig_from_dict(reads)

    @staticmethod
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
                                     qual=read['qual'])
                con = \
                  locate_sequence(seq, location=start, parent=contig)
                contig.consensus = con
            else:
                #reads
                if 'qual' in read:
                    qual = read['qual']
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

