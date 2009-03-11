'''
Created on 2009 mar 9

@author: peio

Script who reads a caf file and takes the information for each
contig. It creates a file for each contig to  easyly use them
'''

from re import match

class CafFile(object):
    ''' This class is used to create and manipulate caf files using an index'''
    def __init__(self, fname):
        
        '''The initialitation.
        keyword arguments:
            fname : caf file name
        '''
        self._fname     = fname
        self._caf_index = {}
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
         
        
        sec_in = False
        type_  = None
        
        fhandler = open(self._fname,'rt')
        rawline  = "Filled"
        while len(rawline) != 0:
            
            prior_tell = fhandler.tell()
            rawline    = fhandler.readline()
            line       = rawline.strip()
            if sec_in:
                if line == "Is_read"  :
                    type_ = "Is_read"
                elif line == "Is_contig": 
                    type_ = "Is_contig"
                
            if match("\w*\s*:\s*\w*", line):
                mode, name = line.split(":")
                mode = mode.strip()
                name = name.strip()
                if mode == "Sequence":
                    sec_in = True
                else:
                    sec_in = False
                
                if name not in self._caf_index.keys():
                    self._caf_index[name] = {}
                
                self._caf_index[name]['name'] = name
                self._caf_index[name]['type'] = type_
                self._caf_index[name][mode]   = prior_tell
                

    def contigs(self):
        '''It returns a generator with the contigs'''
        
        for seq_rec in self._caf_index:
            if self._caf_index[seq_rec]['type'] == 'Is_contig':
                print self._caf_index[seq_rec]['name']
                yield self._get_sec_rec_full(self._caf_index[seq_rec]['name'])

    def reads(self):
        '''It returns a generator with the reads'''
        
        for seq_rec in self._caf_index:
            if self._caf_index[seq_rec]['type'] == 'Is_read':
                yield  self._get_sec_rec_full(self._caf_index[seq_rec]['name'])
  
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
        dna_pos     = self._caf_index[sec_rec_name]['DNA']
        dna_section = self._return_section(dna_pos)
        dna = ''
        for line in dna_section:
            dna += line.strip()
        return dna  

    def _get_base_quality(self, sec_rec_name):
        ''' It returns the base quality array. It needs the sec_rec name '''
         
        base_quality_pos     = self._caf_index[sec_rec_name]['BaseQuality']
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
        return (item[1], item[2], item[3], item[4])
    
    @staticmethod
    def _get_assembled_from(line):
        ''' It reads Assembled_from line and returns 2 elements: The first is
        the read name, and the second one is a tupla with 4 coordinates:
        The order of the tupla is: (contig_start, contig_end, read_start,
         read_end)'''
        item         = line.split(" ")
        return item[1], (item[2], item[3], item[4], item[5])
    
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
        
        sequence_pos     = self._caf_index[sec_rec_name]['Sequence']
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
        
    def _get_sec_rec_full(self, sec_rec_name):
        ''' It returns the complete info of a sec_record. It uses the index 
         to access the file so we do no t need to read the whole file.
         We need the name of the sec record '''
        
        
        
        sec_rec_info = self._get_seq_rec(sec_rec_name) 
        # First we take dna secuence
        sec_rec_info['DNA'] = self._get_dna(sec_rec_name)
        # If BaseQuality is in the sec record       
        if 'BaseQuality' in self._caf_index[sec_rec_name].keys():
            sec_rec_info['BaseQuality'] = self._get_base_quality(sec_rec_name)
        
        return sec_rec_info
    
