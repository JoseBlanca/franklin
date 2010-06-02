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

import itertools, tempfile
from array import array

def write_distribution(fhand, distribution, bin_edges):
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
                        plot_fhand=None, range_=None):
    ''''Given a list of numbers it returns the distribution and it plots the
    histogram'''
    if bins is None:
        bins = 20
    if range_ == (None, None):
        range_ = None
    if labels is None:
        labels = {'title':'histogram', 'xlabel':'values', 'ylabel':'count'}
    #we do the distribution
    result = histogram(numbers, bins=bins, range_=range_)
    if result:
        distrib, bin_edges = result
    else:
        #there is no numbers
        return
    #we write the output
    if distrib_fhand is not None:
        write_distribution(distrib_fhand, distrib, bin_edges)
    #do we have to plot it?
    if plot_fhand is not None:
        draw_histogram(distrib, bin_edges,
                      title=labels['title'], xlabel=labels['xlabel'],
                      ylabel=labels['ylabel'],
                      fhand=plot_fhand)
    return {'distrib':list(distrib), 'bin_edges':list(bin_edges)}

def _range(numbers):
    'Given an iterator with numbers it returns the min and max'
    min_, max_ = None, None
    for number in numbers:
        if min_ is None or min_ > number:
            min_ = number
        if max_ is None or max_ < number:
            max_ = number
    return min_, max_

def histogram(numbers, bins, range_= None):
    '''An alternative implementation to the numpy.histogram.

    The main difference is that it accepts iterators
    '''
    #an iterator for the numbers
    if range_ is None or None in range_ :
        nums, numbers = itertools.tee(numbers)
        calc_range = _range(nums)
        if range_ is None:          #min and max
            range_ = calc_range
        elif range_[0] is None:     #min
            range_[0] = calc_range[0]
        elif range_[1] is None:     #max
            range_[1] = calc_range[1]

    min_, max_ = range_

    if not min_ and not max_:
        #there's no numbers
        return

    #now we can calculate the bin egdes
    distrib_span = max_ - min_
    bin_span     = distrib_span / float(bins)
    bin_edges = [min_ + bin_ * bin_span for bin_ in range(bins + 1)]

    #now we calculate the distribution
    distrib = [0] * bins
    #an iterator for the numbers
    for number in numbers:
        if number > max_ or number < min_:
            continue
        if number == max_:
            #the last value go into the last bin
            bin_ = bins - 1
        else:
            bin_ = int(float(number - min_) / bin_span)
        if bin_ >= 0:
            distrib[bin_] += 1
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


MIN_FREE_MEMORY_PERCENT = 10
MEMORY_CHECK_CYCLES = 10000

def _check_free_memory(percent=None):
    'It checks that we still have more free memory than the given percent'
    if not percent:
        percent = MIN_FREE_MEMORY_PERCENT
    fhand = open('/proc/meminfo', 'r')
    mem_total = int(fhand.readline().split()[1])
    mem_free = int(fhand.readline().split()[1])
    free_percent = mem_free / mem_total * 100
    if percent > free_percent:
        raise RuntimeError('Scarce memory')

class CachedArray(object):
    '''It stores numbers for future use.

    It can work in memory or write it contents to a file.
    '''
    def __init__(self, typecode):
        'The init '
        self._typecode = typecode
        self._array = array(typecode)
        self._cache_fhand = None
        self._max = None
        self._min = None
        self._len = 0
        self._sum = 0
        self._check_in_cycles = MEMORY_CHECK_CYCLES

    def extend(self, items):
        'It adds all items to the store'
        for item in items:
            self.append(item)

    def append(self, item):
        'It appends one item to the store'
        #some statistics
        if self._max is None:
            self._max = item
            self._min = item
        else:
            if self._max < item:
                self._max = item
            if self._min > item:
                self._min = item
        self._sum += item
        self._len += 1

        #are we running low on memory?
        self._check_in_cycles -= 1
        if self._check_in_cycles <= 0 and self._cache_fhand is not None:
            try:
                _check_free_memory()
            except RuntimeError:
                #we are running low on memory, so we store in disk
                self.to_disk()
            self._check_in_cycles = MEMORY_CHECK_CYCLES

        if self._array is not None:
            self._array.append(item)
        else:
            self._cache_fhand.write(str(item) + '\n')

    def _generate_file_items(self):
        'It yields all items from the file cache'
        #let's be sure that all items are in the disk
        if self._cache_fhand:
            self._cache_fhand.flush()

        casts = {'h': int, 'H':int, 'i':int, 'I':int, 'L':int, 'l':int,
                 'f':float }
        cast = casts[self._typecode]
        for line in open(self._cache_fhand.name):
            yield cast(line)

    def _get_items(self):
        'A generator that yields all items'
        if self._array is not None:
            return iter(self._array)
        else:
            return self._generate_file_items()

    def __iter__(self):
        'Part of the iterator protocol'
        return self._get_items()

    def to_disk(self):
        'It saves all items to the disk'
        if self._cache_fhand is not None:
            return  #we are already in disk
        self._cache_fhand = tempfile.NamedTemporaryFile()
        #we store all item in the file
        for item in self._array:
            self._cache_fhand.write(str(item) + '\n')
        self._cache_fhand.flush()
        #we remove the array storage
        del self._array
        self._array = None

    def _get_min(self):
        'It returns the minimum'
        return self._min
    min = property(_get_min)

    def _get_max(self):
        'It returns the maximum'
        return self._max
    max = property(_get_max)

    def _get_sum(self):
        'It returns the sum'
        return self._sum
    sum = property(_get_sum)

    def _get_average(self):
        'It returns the maximum'
        return self._sum / self._len
    average = property(_get_average)

    def __len__(self):
        'It returns the number of items'
        return self._len

    def _get_variance(self):
        'It returns the variance'
        sum_ = 0
        count = 0
        mean = self.average
        for number in self:
            sum_ += (number - mean) ** 2
            count += 1
        if not count:
            return None
        else:
            return (sum_ / float(count))
    variance = property(_get_variance)
