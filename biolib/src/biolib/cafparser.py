'''
Created on 2009 mar 9

@author: peio

Script who reads a caf file and takes the information for each
contig. It creates a file for each contig to  easyly use them
'''

from re import match
from biolib.contig import Contig, locate_sequence
from SeqRecord import SeqRecord

class CafFile(object):
    ''' This class is used to create and manipulate caf files using an index'''
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
        self._caf_file2caf_index()

    def _caf_file2caf_index(self):
        '''It takes a caf file and after reading it it returns an index with 
           section positions. We define section as the paragraph of text where
            each caf data type is represented . example:
                    DNA: name
                    aaaaaaaaaaaaaaa
                    ttttttttttt
                    cccccccccc
         It stores as well if the secuence is a contig or a read
         '''
        qual_index = {}
        seq_index = {}
        dna_index = {}
        type_index = {}
        fhandler = open(self._fname,'rt')
        rawline  = "Filled"
        sec_in = False
        debug_lines = 0
        while len(rawline) != 0:
            
            prior_tell = fhandler.tell()
            rawline    = fhandler.readline()
            line       = rawline.strip()
            
            if debug_lines == 4000000:
                rawline = ''
            debug_lines += 1
            #if match("\w*\s*:\s*\w*", line):
            mode = rawline.split(':', 1)[0].strip().lower()
            if mode in ('dna', 'basequality', 'sequence'):
                mode_, name = line.split(":")
                name = name.strip()
                if mode == "Sequence":
                    sec_in = True
                else:
                    sec_in = False
                if mode == 'dna':
                    dna_index[name] = prior_tell
                elif mode == 'basequality':
                    qual_index[name] = prior_tell
                elif mode == 'sequence':
                    seq_index[name] = prior_tell    

            if sec_in:
                if line == "Is_read":
                    type_index[name] = "Is_read"
                elif line == "Is_contig": 
                    type_index[name] = "Is_contig"
                    
    def contigs(self):
        '''It returns a generator with the contigs'''
        for seq_rec in self._seq_index:
            if self._type_index[seq_rec] == 'Is_contig':
                seq_rec_name = seq_rec
                yield self._get_seq_rec_full(seq_rec_name)

    def reads(self):
        '''It returns a generator with the reads'''
        
        for seq_rec in self._seq_index:
            if self._type_index[seq_rec] == 'Is_read':
                yield  self._get_seq_rec_full(seq_rec)

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

    def _get_dna(self, sec_rec_name):
        ''' It returns the dna secuence in a string.It needs the sec_rec name
        '''
        dna_pos     = self._dna_index[sec_rec_name]
        dna_section = self._return_section(dna_pos)
        dna = ''
        for line in dna_section:
            dna += line.strip()
        return dna
  
    def _get_base_quality(self, sec_rec_name):
        ''' It returns the base quality array. It needs the sec_rec name '''
        base_quality_pos     = self._qual_index[sec_rec_name]
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
        the read name, and the second one is a tupla with 4 coordinates:
        The order of the tupla is: (contig_start, contig_end, read_start,
         read_end)'''
        item         = line.split(" ")
        return item[1], (int(item[2]), int(item[3]), int(item[4]), int(item[5]))

    def _get_seq_rec(self, sec_rec_name):
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
                if read_name not in reads.keys():
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
        
        seq_rec_info = self._get_seq_rec(seq_rec_name) 
        # First we take dna secuence
        seq_rec_info['DNA'] = self._get_dna(seq_rec_name)
        # If BaseQuality is in the sec record
        if seq_rec_name in self._qual_index:
            seq_rec_info['BaseQuality'] = self._get_base_quality(seq_rec_name)
        
        return seq_rec_info
    
    def read_record2read(self, seq_rec_name):
        ''' This wraper takes the dictionary taked from each sec record
        and put them into a biopython seq_rec objetc'''
        
        if self._type_index[seq_rec_name] == 'Is_read':
            
            seq_info    = self._get_seq_rec(seq_rec_name)
            seq_dna     = self._get_dna(seq_rec_name)
            seq_quality = self._get_base_quality(seq_rec_name)
     
            seq_rec  = SeqRecord(seq = seq_dna, id = seq_rec_name, \
                                 name = seq_rec_name )
            seq_rec.letter_annotations['quality'] = seq_quality
            seq_rec.annotations =  seq_info
            return seq_rec
        else:
            print seq_rec_name
            raise RuntimeError ('This sec record it supposed to be a read\
                                 record and is not')
    @staticmethod  
    def _correct_minus(reads):
        ''' It corrects the problem of the minus(-)  coordenates. This function
        returns a int that is de maximun minus numer.example:
                   -2101234567890
            contig    aaaaaaaa
            read1   cccccc
            read2      eeeeee
            In this case the maximun minus number is 2 and we use it to 
            move the other seqs
            ''' 
        correction = 0
        for read in reads:
            contig_start = reads[read][0][0]
            read_start   = reads[read][0][2]
            diff         = contig_start - read_start
            if diff < correction:
                correction = diff
        return abs(correction)
                
    def contig_record2contig(self, contig_name):
        ''' This wraper takes the dictionary taked from each sec record
         and put them into a contig objetc or a sec_record object'''

        contig_info = self._get_seq_rec(contig_name)
        contig_dna  = self._get_dna(contig_name)
        if contig_info['type'] == 'Is_contig':
            reads      = contig_info['reads']
            correction = self._correct_minus(reads)
           
            if len(contig_dna) == 0:
                contig     = Contig()
            else:
                consensus = SeqRecord(seq=contig_dna, id='Consensus',
                                      name = contig_name)
                consensus  = locate_sequence(sequence = consensus, \
                                             location = correction)
                contig     = Contig(consensus=consensus)
                
            for read in reads:
                read_sections = reads[read]
                if len(read_sections) > 1:
                    #We it's unpadded there is more than one list of coordentes 
                    raise RuntimeError('This is an unppadded sec record and \
                     we do not admit them. Use caf tools to convert it to \
                     padded')
                else:
                    contig_start  = int(read_sections[0][0]) + correction
                    read_start    = int(read_sections[0][2]) 
                    contig_start -= read_start

                #Data to fill the contig
                seq_rec = self.read_record2read(read)
                sec_strand = seq_rec.annotations['Strand']
#                print seq_rec.annotations['Clipping']
                mask = seq_rec.annotations['Clipping'].split(" ")[1:]
                
                if sec_strand.lower() == 'forward':
                    strand  = 1
                    forward = True
                else:
                    strand = -1
                    forward = False
                
                contig.append_to_location(sequence=seq_rec, \
                                          start=contig_start,\
                                          strand=strand, forward=forward, \
                                          mask=mask)
        return contig
                 
    
            
        
         
         
        