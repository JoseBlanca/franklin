'''Some functions to do statistics on sequences.'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

import numpy
from franklin.utils.collections_ import FileCachedList
import itertools

PLOT_LABELS = {'masked_seq_distrib' :{
                                 'title': 'Masked sequence length distribution',
                                 'xlabel':'Masked length',
                                 'ylabel':'Number of sequences'},
                'seq_length_distrib':{'title':'Sequence length distribution',
                                      'xlabel':'Sequence length',
                                      'ylabel': 'Number of sequences'},
                'qual_distrib'      :{'title':'Sequence quality distribution',
                                      'xlabel':'quality',
                                      'ylabel':'Number of bp'},
                'contig_read_distrib':{'title':'Contig reads distribution',
                                      'xlabel':'number of reads',
                                      'ylabel':'number of contigs'},
                'contig_coverage_distrib':{
                                      'title':'Contig coverage distribution',
                                      'xlabel':'coverage',
                                      'ylabel':'number of positions'},
                }

def _write_distribution(fhand, distribution, bin_edges):
    '''It writes the given distribution in a file.

    It requires a list with the distribution and a list with the bin_edges.
    len(bin_edges) == len(distribution) + 1
    '''
    for name, numbers in (('values', distribution), ('bin_edges', bin_edges)):
        fhand.write(name + ': ')
        fhand.write(','.join((str(item) for item in numbers)))
        fhand.write('\n')
    fhand.flush()

def create_distribution(numbers, labels=None, distrib_fhand=None, bins=None,
                        plot_fhand=None, range_=None, low_memory=True):
    ''''Given a list of numbers it returns the distribution and it plots the
    histogram'''
    if bins is None:
        bins = 20
    if range_ == (None, None):
        range_ = None
    if labels is None:
        labels = {'title':'histogram', 'xlabel':'values', 'ylabel':'count'}
    #we do the distribution
    if low_memory:
        distrib, bin_edges = histogram(numbers, bins=bins, range_=range_)
    else:
        distrib, bin_edges = numpy.histogram(numbers, bins=bins, new=True,
                                         range=range_)
    #we write the output
    if distrib_fhand is not None:
        _write_distribution(distrib_fhand, distrib, bin_edges)
    #do we have to plot it?
    if plot_fhand is not None:
        draw_histogram(distrib, bin_edges,
                     title=labels['title'], xlabel=labels['xlabel'],
                     ylabel=labels['ylabel'],
                     fhand=plot_fhand)
    return {'distrib':list(distrib), 'bin_edges':list(bin_edges)}

def _masked_sequence_lengths(sequences, low_memory):
    'It returns a list with the sequences lengths'
    if low_memory:
        lengths = FileCachedList(float)
    else:
        lengths = []
    for seq in sequences:
        length = 0
        for letter in str(seq.seq):
            if letter.islower():
                length += 1
        lengths.append(length)
    return lengths

def _sequence_lengths(sequences, low_memory):
    'It returns a list with the sequences lengths'
    if low_memory:
        lengths = FileCachedList(float)
    else:
        lengths = []
    for seq in sequences:
        lengths.append(len(seq))
    return lengths

def _sequence_qualitities(sequences, low_memory):
    'It returns a list with the sequence qualities'
    if low_memory:
        qualities = FileCachedList(float)
    else:
        qualities = []
    for seq in sequences:
        # This statistic can be run when there is no quality. In this cases it
        # should exit without warnings
        if seq.qual is None:
            continue
        for value in seq.qual:
            qualities.append(value)
    return qualities

def _contig_reads(contigs, low_memory):
    'It returns a list with the number of reads per contig'
    if low_memory:
        values = FileCachedList(float)
    else:
        values = []
    for contig in contigs:
        values.append(len(contig))
    return values

def _contig_to_read_list(contig):
    ''' It takes a contig class object and it fill a list with the reads.
    All the reads are '''
    reads = []
    for read in contig:
        reads.append(str(read))
    return reads

def _contig_coverages(contigs, low_memory):
    'It returns a list with the number of reads per position in the contigs'
    if low_memory:
        values = FileCachedList(float)
    else:
        values = []
    for contig in contigs:
        contig = _contig_to_read_list(contig)
        ncols  = len(contig[0])
        for col in range(ncols):
            coverage = 0
            for read in contig:
                try:
                    base = read[col]
                    if base and not base.isspace():
                        coverage += 1
                    #pylint: disable-msg=W0704
                except IndexError:
                    pass
            values.append(coverage)
    return values

#pylint:disable-msg=W0613
def seq_distrib(kind, sequences, distrib_fhand=None, plot_fhand=None,
                range_=None, low_memory=False):
    'It returns the stat function depending on the stat it calculates(type)'
    stats = {'masked_seq_distrib'     : _masked_sequence_lengths,
             'seq_length_distrib'     : _sequence_lengths,
             'qual_distrib'           : _sequence_qualitities,
             'contig_read_distrib'    : _contig_reads,
             'contig_coverage_distrib': _contig_coverages,
             }
    values = stats[kind](sequences, low_memory)
    labels = PLOT_LABELS[kind]
    return create_distribution(values, labels, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand, range_=range_,
                                low_memory=low_memory)

def seq_distrib_diff(seqs1, seqs2, kind, distrib_fhand=None, plot_fhand=None,
                     low_memory=True):
    'It return the difference between different distributions of the given type'

    #get labesl depending on the stat type
    labels = PLOT_LABELS[kind]

    #the values for boths sequences
    values_functs = {
             'masked_seq_distrib': _masked_sequence_lengths,
             'seq_length_distrib': _sequence_lengths,
             'qual_distrib'      : _sequence_qualitities}
    values1 = values_functs[kind](seqs1, low_memory=low_memory)
    values2 = values_functs[kind](seqs2, low_memory=low_memory)

    #the range
    if low_memory:
        num_iters1, num_iters2 = values1.items(), values2.items()
    else:
        num_iters1, num_iters2 = iter(values1), iter(values2)

    min_, max_ = None, None
    for values in (num_iters1, num_iters2):
        for number in values:
            if min_ is None or min_ > number:
                min_ = number
            if max_ is None or max_ < number:
                max_ = number
    range_ = (min_, max_)

    #Now I calculate the distibution with the same range, in order to be able
    # to get the difference
    if low_memory:
        num_iters1, num_iters2 = values1.items(), values2.items()
    else:
        num_iters1, num_iters2 = iter(values1), iter(values2)

    distrib1 = create_distribution(num_iters1,  range_=range_,
                                    low_memory=low_memory)
    distrib2 = create_distribution(num_iters2,  range_=range_,
                                    low_memory=low_memory)

    diff_distrib   = []
    diff_bin_edges = distrib1['bin_edges']
    # now a subtract distrib1 from distrib2
    for i in range(len(distrib1['distrib'])):
        diff = distrib1['distrib'][i] - distrib2['distrib'][i]
        diff_distrib.append(diff)

    #we write the output
    if distrib_fhand is not None:
        _write_distribution(distrib_fhand, diff_distrib, diff_bin_edges)
    #do we have to plot it?
    if plot_fhand is not None:
        draw_histogram(diff_distrib, diff_bin_edges,
                     title=labels['title'], xlabel=labels['xlabel'],
                     ylabel=labels['ylabel'],
                     fhand=plot_fhand)
    return {'distrib':diff_distrib, 'bin_edges':diff_bin_edges}

def _range(numbers):
    'Given an iterator with numbers it returns the min and max'
    min_, max_ = None, None
    for number in numbers:
        if min_ is None or min_ > number:
            min_ = number
        if max_ is None or max_ < number:
            max_ = number
    return min_, max_

def general_seq_statistics(sequences, low_memory=True):
    '''Given a sequence iterator it calculates some general statistics.

    The statistics will be written into the given distrib_fhand (if given) and
    a dict with them will be returned.
    It calculates the total sequence length, the average sequence length, the
    total masked sequence length, and the number of sequences.
    '''
    has_qual = False
    if low_memory:
        lengths        = FileCachedList(int)
        masked_lengths = FileCachedList(int)
        qualities      = FileCachedList(int)
    else:
        lengths        = []
        masked_lengths = []
        qualities      = []

    n_seqs = 0
    total_len = 0
    total_masked_len = 0
    for seq in sequences:
        n_seqs += 1
        length = len(seq)
        lengths.append(length)
        total_len += length
        #low_memory False, is just one sequence
        masked_length = _masked_sequence_lengths([seq], low_memory=False)[0]
        total_masked_len += masked_length
        masked_lengths.append(masked_length)
        if seq.qual:
            has_qual = True
            qualities.extend(seq.qual)

    stats = {}
    stats['num_sequences']     = n_seqs
    stats['seq_length']        = total_len
    stats['masked_seq_length'] = total_masked_len
    mean_quality     = None
    quality_variance = None
    min_quality      = None
    max_quality      = None
    if low_memory:
        stats['seq_length_average']         = _average(lengths.items())
        stats['seq_length_variance']        = _variance(lengths.items(),
                                                 stats['seq_length_average'])
        stats['min_seq_length'], stats['max_seq_length'] = \
                                                         _range(lengths.items())
        stats['masked_seq_length_average']  = _average(masked_lengths.items())
        stats['masked_seq_length_variance'] = _variance(masked_lengths.items(),
                                            stats['masked_seq_length_average'])
        if has_qual:
            mean_quality     = _average(qualities.items())
            quality_variance = _variance(qualities.items(), mean_quality)
            min_quality, max_quality = _range(qualities.items())
    else:
        stats['seq_length_average']         = numpy.average(lengths)
        stats['seq_length_variance']        = numpy.var(lengths)
        stats['min_seq_length'], stats['max_seq_length'] = _range(lengths)
        stats['masked_seq_length_average']  = numpy.average(masked_lengths)
        stats['masked_seq_length_variance'] = numpy.var(masked_lengths)
        if has_qual:
            mean_quality     = _average(qualities)
            quality_variance = _variance(qualities, mean_quality)
            min_quality, max_quality = _range(qualities)
    stats['mean_quality']     = mean_quality
    stats['quality_variance'] = quality_variance
    stats['min_quality']      = min_quality
    stats['max_quality']      = max_quality
    return stats

def _average(numbers):
    'Given an iterator with numbers it calculates the average'
    sum_ = 0
    count = 0
    for number in numbers:
        sum_ += number
        count += 1
    return sum_ / float(count)

def _variance(numbers, mean):
    'Given an iterator with numbers it calculates the variance'
    sum_ = 0
    count = 0
    for number in numbers:
        sum_ += (number - mean) ** 2
        count += 1
    return (sum_ / float(count))

def general_contig_statistics(contigs, distrib_fhand=None, low_memory=True):
    'It calculates some contigs statistics'
    if low_memory:
        coverages = FileCachedList(float)
        reads     = FileCachedList(float)
    else:
        coverages = []
        reads     = []

    n_contigs = 0
    for contig in contigs:
        n_contigs += 1
        #low_memory False becasue is just one contig
        coverages.extend(_contig_coverages([contig], low_memory=False))
        reads.extend(_contig_reads([contig], low_memory=False))

    stats = {}
    stats['number_contigs'] = n_contigs
    if low_memory:
        stats['mean_number_reads']     = _average(reads.items())
        stats['number_reads_variance'] = _variance(reads.items(),
                                                   stats['mean_number_reads'])
        stats['mean_coverage']         = _average(coverages.items())
        stats['coverage_variance']     = _variance(coverages.items(),
                                                   stats['mean_coverage'])
    else:
        stats['mean_number_reads']     = numpy.average(reads)
        stats['number_reads_variance'] = numpy.var(reads)
        stats['mean_coverage']         = numpy.average(coverages)
        stats['coverage_variance']     = numpy.var(coverages)

    return stats

def histogram(numbers, bins, range_=None):
    '''An alternative implementation to the numpy.histogram.

    The main difference is that this implementation can use a CachedFileList to
    save memory
    '''
    #pylint: disable-msg=W0622
    #is this a FileCachedList?
    cached = False
    if 'items' in dir(numbers):
        cached = True
    #an iterator for the numbers
    if cached:
        num_iter  = numbers.items()
        num_iter2 = numbers.items()
    else:
        num_iter, num_iter2 = itertools.tee(iter(numbers), 2)
    if range_ is None:
        min_, max_ = _range(num_iter)
    else:
        min_, max_ = range_[0], range_[1]
    if min_ is None:
        min_ = _range(num_iter)[0]
    if max_ is None:
        max_ = _range(num_iter)[1]

    #now we can calculate the bin egdes
    distrib_span = max_ - min_
    bin_span     = distrib_span / float(bins)
    bin_edges = [min_ + bin * bin_span for bin in range(bins + 1)]

    #now we calculate the distribution
    distrib = [0] * bins
    #an iterator for the numbers
    for number in num_iter2:
        if number > max_ or number < min_:
            continue
        if number == max_:
            #the last value go into the last bin
            bin = bins - 1
        else:
            bin = int(float(number - min_) / bin_span)
        if bin >= 0:
            distrib[bin] += 1
    return (distrib, bin_edges)

IMPORTED_MATPLOTLIB = None

def draw_histogram(values, bin_edges, title=None, xlabel= None, ylabel=None,
                   fhand=None):
    'It draws an histogram and if the fhand is given it saves it'
    global IMPORTED_MATPLOTLIB
    if IMPORTED_MATPLOTLIB is None:
        import matplotlib
        matplotlib.use('AGG')
        IMPORTED_MATPLOTLIB = matplotlib

    import matplotlib.pyplot as plt

    fig = plt.figure()
    axes = fig.add_subplot(111)
    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)
    if title:
        axes.set_title(title)

    xvalues = range(len(values))

    axes.bar(xvalues, values)

    #the x axis label
    xticks_pos = [value + 0.5 for value in xvalues]

    left_val = None
    right_val = None
    xticks_labels = []
    for value in bin_edges:
        right_val = value
        if left_val:
            xticks_label = (left_val + right_val) / 2.0
            if xticks_label >= 10:
                fmt = '%d'
            elif xticks_label >= 0.1 and xticks_label < 10:
                fmt = '%.1f'
            elif xticks_label < 0.1:
                fmt = '%.1e'
            xticks_label = fmt % xticks_label
            xticks_labels.append(xticks_label)
        left_val = right_val

    #we don't want to clutter the plot
    xticks_pos = xticks_pos[::2]
    xticks_labels = xticks_labels[::2]
    plt.xticks(xticks_pos, xticks_labels)

    if fhand is None:
        plt.show()
    else:
        plt.savefig(fhand)

def calculate_read_coverage(pileup, distrib_fhand=None, plot_fhand=None,
                            range_=None):
    '''Given a sam pileup file it returns the coverage distribution.

    The coverage shows how many times the bases has been read.
    '''
    coverages = FileCachedList(int)
    for line in pileup:
        if line.isspace():
            continue
        position_cov = line.split()[3]
        coverages.append(position_cov)
    #now the distribution
    return create_distribution(coverages,
                               labels={'title':'Read coverage distribution',
                                      'xlabel':'coverage',
                                      'ylabel': 'Number of positions'},
                               distrib_fhand=distrib_fhand,
                               plot_fhand=plot_fhand,
                               range_=range_, low_memory=True)

def _color_by_index(index, kind='str'):
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
    #in those occasions
    global IMPORTED_MATPLOTLIB
    if IMPORTED_MATPLOTLIB is None:
        import matplotlib
        matplotlib.use('AGG')
        IMPORTED_MATPLOTLIB = matplotlib

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
        scatters[scat_index]['color'] = _color_by_index(color_index,
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

