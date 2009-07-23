'''
Created on 2009 api 30

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


import tempfile, shutil
from uuid import uuid4
import os
import math

from biolib.seqs import SeqWithQuality

from Bio import  SeqIO

def float_lists_are_equal(list1, list2):
    'Given two lists it checks that all floats are equal'
    for num1, num2 in zip(list1, list2):
        assert floats_are_equal(num1, num2)

def floats_are_equal(num1, num2):
    'Given two numbers it returns True if they are similar'
    if num1 == 0.0:
        if num2 == 0.0:
            return True
        else:
            return False
    log1 = math.log(float(num1))
    log2 = math.log(float(num2))
    return abs(log1 - log2) < 0.01

class NamedTemporaryDir(object):
    '''This class creates temporary directories '''
    def __init__(self):
        '''It initiates the class.'''
        self._name = tempfile.mkdtemp()

    def get_name(self):
        'Returns path to the dict'
        return self._name
    name = property(get_name)
    def close(self):
        '''It removes the temp dir'''
        if os.path.exists(self._name):
            shutil.rmtree(self._name)

    def __del__(self):
        '''It removes de temp dir when instance is removed and the garbaje
        colector decides it'''
        self.close()

def fasta_str(seq, name):
    'Given a sequence it returns a string with the fasta'
    fasta_str_ = ['>']
    fasta_str_.append(name)
    fasta_str_.append('\n')
    fasta_str_.append(str(seq))
    fasta_str_.append('\n')
    return ''.join(fasta_str_)

def _get_seq_name(seq):
    'Given a sequence and its default name it returns its name'
    try:
        name = seq.name
    except AttributeError:
        name = None

    if name is None:
        name  = str(uuid4())
    return name

def temp_fasta_file(seqs, write_qual=False):
    '''Given a Seq and its default name it returns a fasta file in a
    temporary file. If the seq is a SeqWithQuality you can ask a qual fasta
    file'''


    try:
        # Is seqs an seq or an iter??
        #pylint:disable-msg=W0104
        seqs.name
        seqs = [seqs]
        #pylint:disable-msg=W0704
    except AttributeError:
        pass

    fhand_seq = tempfile.NamedTemporaryFile(suffix='.fasta')
    if write_qual:
        fhand_qual = tempfile.NamedTemporaryFile(suffix='.fasta')

    for seq in seqs:
        name = _get_seq_name(seq)

        if write_qual:
            try:
                quality = seq.qual
            except AttributeError:
                msg = 'Sequence must be a SeqWithQuality instance'
                raise AttributeError(msg)
            if quality is not None:
                quality = [str(qual) for qual in quality]
                fhand_qual.write(fasta_str(' '.join(quality), name))
            else:
                raise AttributeError('Quality can not be empty')
        fhand_seq.write(fasta_str(seq, name))

    fhand_seq.flush()
    fhand_seq.seek(0)
    if write_qual:
        fhand_qual.flush()
        fhand_qual.seek(0)
        return fhand_seq, fhand_qual
    else:
        return fhand_seq

def temp_multi_fasta_file(seqs):
    '''It creates a temporari fasta file using a list of sequences'''
    fileh = tempfile.NamedTemporaryFile(suffix='.fasta')
    for seq in seqs:
        name = _get_seq_name(seq)
        fileh.write(fasta_str(seq, name))
    fileh.flush()
    fileh.seek(0)
    return fileh

def create_temp_fasta_files(seq1, seq2):
    'It returns two temporal fasta files.'
    #we create two temp files
    fileh1 = temp_fasta_file(seq1)
    fileh2 = temp_fasta_file(seq2)
    return fileh1, fileh2

def get_start_end(location):
    '''It accepts an int, Location or tuple and it returns the start, end,
    forward and strand.'''
    #int
    if isinstance(location, int):
        start = location
        end = location
    #tuple
    elif isinstance(location, tuple):
        start = location[0]
        end = location[1]
    #location
    else:
        start = location.start
        end = location.end
    return start, end

def parse_fasta(seq_fhand, qual_fhand=None):
    '''It returns the fasta file content giving a file hanler'''
    seq  = []
    qual = None
    name, description, seq = get_content_from_fasta(seq_fhand)
    if qual_fhand is not None:
        #pylint: disable-msg=W0612
        name_, descrip, qual = get_content_from_fasta(qual_fhand, kind='qual')
    if not seq and not qual:
        return None
    else:
        return SeqWithQuality(seq=seq, qual=qual, name=name,
                          description=description)

def get_content_from_fasta(fhand, kind='seq'):
    '''It returns the seq/qual from a fasta file, it need a fhand'''
    seq         = []
    name        = None
    description = None
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            items = line.split(' ', 1)
            name  = items[0][1:]
            try:
                description = items[1]
            except IndexError:
                description = None
        else:
            seq.append(line)
    if kind == 'seq':
        seq = ''.join(seq)
    else:
        qual = []
        for item in seq:
            qual.extend(item.split())
        qual = [int(qitem) for qitem in qual]
        seq = qual
    return name, description, seq

def remove_from_orf(orf_dna, orf_prot, aminoacid='X'):
    ''' It removes an aminoaacid from dna and protein seq'''
    dna  = []
    prot = []
    pos       = 0
    for letter in orf_prot:
        if letter.upper() != aminoacid:
            prot.append(letter)
            dna.append(orf_dna[pos:pos + 3])
        pos += 3
    return "".join(dna), "".join(prot)

def _remove_atributes_to_tag(tag):
    '''It removees atributes to a xml tag '''
    mod_tag = "".join(tag).split(' ')
    if len(mod_tag) >1:
        return mod_tag[0] + '>'
    else:
        return mod_tag[0]

def _get_xml_header(fhand, tag):
    '''It takes the header of the xml file '''
    fhand.seek(0, 2)
    end_file = fhand.tell()

    fhand.seek(0, 0)
    header      = []
    current_tag = []
    listed_tag = '<' + tag + '>'
    while True:
        if end_file <= fhand.tell():
            raise ValueError('End Of File. Tag Not found')

        letter = fhand.read(1)
        if letter == '<':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        else:
            header.append(letter)
        if letter == '>':
            mod_tag = _remove_atributes_to_tag(current_tag)
            if listed_tag == mod_tag:
                return  ''.join(header)
            else:
                header.extend(current_tag)
                current_tag = []
            in_tag = False
def _get_xml_tail(fhand, tag):
    '''It takes the tail of the xml file '''
    in_tag = False
    tail = []
    current_tag  = []
    fhand.seek(-1, 2)
    listed_tag = list('</'+ tag +'>')
    listed_tag.reverse()
    while True:
        if fhand.tell() == 0:
            raise ValueError('Start Of File. Tag Not found')
        letter = fhand.read(1)
        if letter == '>':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        else:
            tail.append(letter)

        if letter == '<':
            if current_tag == listed_tag:
                tail.reverse()
                return "".join(tail)
            else:
                tail.extend(current_tag)
                current_tag = []
            in_tag = False

        fhand.seek(-2, 1)

def xml_itemize(fhand, tag):
    '''It takes a xml file and it chunks it by the given key. It adds header if
    exists to each of the pieces. It is a generator'''
    fhand.seek(0, 2)
    end_file = fhand.tell()

    header = _get_xml_header(fhand, tag)
    tail   = _get_xml_tail(fhand, tag)
    section      = []
    current_tag  = []
    listed_tag_s = '<' + tag + '>'
    listed_tag_e = '</' + tag + '>'
    in_tag, in_section = False, False

    fhand.seek(0, 0)

    while True:
        if end_file <= fhand.tell():
            break
        letter = fhand.read(1)
        if letter == '<':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        if in_section:
            section.append(letter)
        if letter == '>':
            if listed_tag_s == _remove_atributes_to_tag(current_tag):
                in_section = True
                section.extend(current_tag)
            elif listed_tag_e == _remove_atributes_to_tag(current_tag):
                yield  header + "".join(section) + tail
                section = []
                in_section = False
            in_tag = False
            current_tag = []

#xmap and xfilter a taken from
#http://code.activestate.com/recipes/66448/

def color_by_index(index, kind='str'):
    'Given an int index it returns a color'
    colors = [{'black'        :(0x00, 0x00, 0x00)},
              {'green'        :(0x00, 0x80, 0x00)},
              {'deep_sky_blue':(0x00, 0xbf, 0xff)},
              {'indigo'       :(0x4b, 0x00, 0x82)},
              {'maroon'       :(0x80, 0x00, 0x00)},
              {'blue_violet'  :(0x8a, 0x2b, 0xe2)},
              {'pale_green'   :(0x98, 0xfb, 0x98)},
              {'sienna'       :(0xa0, 0x52, 0x22)},
              {'medium_orchid':(0xba, 0x55, 0xd3)},
              {'rosy_brown'   :(0xbc, 0x8f, 0x8f)},
              {'chocolate'    :(0xd2, 0x69, 0x1e)},
              {'crimson'      :(0xdc, 0x14, 0x3c)},
              {'dark_salmon'  :(0xe9, 0x96, 0x7a)},
              {'khaki'        :(0xf0, 0xe6, 0x8c)},
              {'red'          :(0xff, 0x00, 0x00)},
              {'blue'         :(0x00, 0x00, 0xff)},
              {'lime'         :(0x00, 0xff, 0x00)},
             ]
    color = colors[index].values()[0]
    #rgb str
    if kind == 'str':
        color = ' '.join([str(channel) for channel in list(color)])
    elif kind == 'rgb_float':
        color = [float(channel)/255.0 for channel in list(color)]
    return color

def draw_scatter(x_axe, y_axe, names=None, groups_for_color=None,
                 groups_for_shape=None, title=None, xlabel= None,
                 ylabel=None, fhand=None):
    '''It draws an scatter plot.

    x_axe and y_axe should be two lists of numbers. The names should be a list
    of the same lenght and will be applied to every data point drawn. The
    groups_for_shape and color should be list of the same length and will be
    used to calculate the color. Every point that belongs to the same group
    will be drawn with the same color.
    If an fhand is given the scatter plot will be saved.
    '''
    #in some circunstances matplot lib could generate this error
    #Failed to create %s/.matplotlib; consider setting MPLCONFIGDIR to a
    #writable directory for matplotlib configuration data
    #in that case we don't know how to use matplotlib, it would require
    #to set the MPLCONFIGDIR variable, but we can't do that in the
    #current shell, so the matplotlib greatness wouldn't be available
    #in those ocasions
    import matplotlib
    matplotlib.use('AGG')
    import matplotlib.pyplot as plt


    fig = plt.figure()
    axes = fig.add_subplot(111)
    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)
    if title:
        axes.set_title(title)
    #text labels
    if names is not None:
        max_x = max(x_axe)
        min_x = min(x_axe)
        x_text_offset = float(max_x - min_x) / 40.0
        for name, x_value, y_value in zip(names, x_axe, y_axe):
            axes.text(x_value + x_text_offset, y_value, name)

    if names is None:
        #all belong to the same group
        names = (None,) * len(x_axe)
    if groups_for_color is None:
        #all belong to the same group
        groups_for_color = (None,) * len(x_axe)
    if groups_for_shape is None:
        #all belong to the same group
        groups_for_shape = (None,) * len(x_axe)
    #now I want the x, y values divided by color and shape
    scatters = {}
    scat_indexes = []
    for x_value, y_value, name, color, shape in zip(x_axe, y_axe, names,
                                                    groups_for_color,
                                                    groups_for_shape):
        scat_index = str(color) + str(shape)
        if scat_index not in scatters:
            scat_indexes.append(scat_index)
            scatters[scat_index] = {}
            scatters[scat_index]['x'] = []
            scatters[scat_index]['y'] = []
            scatters[scat_index]['names'] = []
            scatters[scat_index]['color'] = color
            scatters[scat_index]['shape'] = shape
        scatters[scat_index]['x'].append(x_value)
        scatters[scat_index]['y'].append(y_value)
        scatters[scat_index]['names'].append(name)
    #which color every scatter should use?
    colors = []
    shapes = []
    avail_shapes = ['o', 's', '^', 'd', 'p', '+', 'x', '<', '>', 'v', 'h']
    for scat_index in scat_indexes:
        color = scatters[scat_index]['color']
        shape = scatters[scat_index]['shape']
        if color not in colors:
            colors.append(color)
        if shape not in shapes:
            shapes.append(shape)
        color_index = colors.index(color)
        shape_index = shapes.index(shape)
        scatters[scat_index]['color'] = color_by_index(color_index,
                                                       kind='rgb_float')
        scatters[scat_index]['shape'] = avail_shapes[shape_index]

    #now the drawing
    for scat_index in scatters:
        scat = scatters[scat_index]
        axes.scatter(scat['x'], scat['y'], c=scat['color'],
                     marker=scat['shape'], s=60)
    if fhand is None:
        plt.show()
    else:
        fig.savefig(fhand)

def guess_seq_file_format(fhand):
    'Given a sequence file it returns its format'
    fhand.seek(0)
    line = fhand.readline()
    if line[0] == '>':
        item = fhand.readline().strip()[0]
        if item.isdigit():
            format_ = 'qual'
        else:
            format_ = 'fasta'
    elif line.split()[0] == 'LOCUS':
        format_ = 'genbank'

    else:
        raise ValueError('Unknown sequence file format for : ' + fhand.name)
    fhand.seek(0)
    return format_

def seqs_in_file(seq_fhand, qual_fhand=None):
    'It yields a seqrecord for each of the secuences found in the seq file'
    seq_file_format = guess_seq_file_format(seq_fhand)

    seq_iter    = SeqIO.parse(seq_fhand, seq_file_format)
    if qual_fhand is None:
        qual_iter = None
    else:
        qual_file_format = guess_seq_file_format(qual_fhand)
        qual_iter = SeqIO.parse(qual_fhand, qual_file_format)

    for seqrec in seq_iter:
        if qual_iter is None:
            qual = None
            qual_name = None
        else:
            qual_sec_record = qual_iter.next()
            qual = qual_sec_record.letter_annotations["phred_quality"]
            qual_name = qual_sec_record.name
        seq         = seqrec.seq
        name        = seqrec.name

        if qual_name and qual_name != name:
            msg = 'Seqs and quals not in the same order: %s, %s' % (name ,
                                                                    qual_name)
            raise RuntimeError(msg)
        description = seqrec.description
        annotations = seqrec.annotations
        yield SeqWithQuality(seq=seq, qual=qual, name=name,
                            description=description, annotations=annotations)

def get_safe_fname(directory, prefix, suffix):
    '''It looks if the name exits in this directory and adds a number to the end
    to be a threath same name'''

    if not os.path.exists('%s/%s.%s' % (directory, prefix, suffix)):
        return '%s/%s.%s' % (directory, prefix, suffix)
    for i in range(10000):
        if not os.path.exists('%s/%s.%d.%s' %(directory, prefix, i,
                                                     suffix)):
            return '%s/%s.%d.%s' % (directory, prefix, i, suffix)

    msg  = 'There are more than 100000 outputs of this kind is this directory'
    raise ValueError(msg)





