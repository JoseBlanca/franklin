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

from __future__  import division

import itertools, tempfile, sys, random
from array import array

try:
    import numpy
except ImportError:
    pass

SAMPLE_LENGTH = 10000

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
                        plot_fhand=None, range_=None, summary_fhand=None,
                        calculate_freqs=False):
    ''''Given a list of numbers it returns the distribution and it plots the
    histogram'''
    if bins is None:
        bins = 20
    if range_ == (None, None):
        range_ = None
    default_labels = {'title':'histogram', 'xlabel':'values', 'ylabel':'count',
                  'minimum':'minimum', 'maximum':'maximum', 'average':'average',
                  'variance':'variance', 'sum':'sum', 'items':'items'}
    if labels is None:
        labels = default_labels
    else:
        for label, value in default_labels.items():
            if label not in labels:
                labels[label] = value
    #we do the distribution
    result = histogram(numbers, bins=bins, range_=range_,
                       calculate_freqs=calculate_freqs)
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
    #now we write some basic stats
    format_num = lambda x: str(x) if isinstance(x, int) else '%.2f' % x
    if summary_fhand:
        msg = 'Statistics for %s\n' % labels['title']
        summary_fhand.write(msg)
        summary_fhand.write('-' * len(msg) + '\n')
        summary_fhand.write('%s: %s\n' % (labels['minimum'],
                                          format_num(numbers.min)))
        summary_fhand.write('%s: %s\n' % (labels['maximum'],
                                          format_num(numbers.max)))
        summary_fhand.write('%s: %s\n' % (labels['average'],
                                          format_num(numbers.average)))
        summary_fhand.write('%s: %s\n' % (labels['variance'],
                                          format_num(numbers.variance)))
        if labels['sum'] is not None:
            summary_fhand.write('%s: %s\n' % (labels['sum'],
                                              format_num(numbers.sum)))
        summary_fhand.write('%s: %s\n' % (labels['items'], len(numbers)))

        summary_fhand.write('\n')
    return {'distrib':list(distrib), 'bin_edges':list(bin_edges)}

def _range(numbers):
    'Given an iterator with numbers it returns the min and max'
    min_, max_ = None, None
    for number in numbers:
        if min_ is None or min_ > number:
            min_ = number
        if max_ is None or max_ < number:
            max_ = number
    return [min_, max_]

def _calculate_range_for_histogram(numbers, range_, remove_outliers):
    'Given a list of numbers it returns a reasonable range'

    if range_ is not None:
        min_, max_ = range_
    else:
        min_, max_ = None, None
    if min_ is not None and max_ is not None:
        return min_, max_

    #min or max is not specified, so we have to calculate them
    if 'min' in dir(numbers):   #is a Storage
        real_range = [numbers.min, numbers.max]
    else:
        nums, numbers = itertools.tee(numbers)
        real_range = _range(nums)

    if remove_outliers:
        percents = [5, 95]
        if 'min' in dir(numbers):   #is a Storage
            outliers_range = numbers.get_sample_percentiles(percents)
        else:
            #we take a sample
            if '__len__' not in dir(numbers):
                #if we have an iterator we have to tee
                nums, numbers = itertools.tee(numbers)
                sample = []
                for index in range(SAMPLE_LENGTH):
                    try:
                        sample.append(nums.next())
                    except StopIteration:
                        break
            else:
                sample = numbers[:SAMPLE_LENGTH]
            outliers_range = _calculate_percentiles(sample, percents)
        #if the limits are close to the real limits we use the real limits
        closeness = outliers_range[1] - outliers_range[0] / 10
        if outliers_range[0] - closeness > real_range[0]:
            real_range[0] = outliers_range[0]
        if outliers_range[1] + closeness < real_range[1]:
            real_range[1] = outliers_range[1]

    if min_ is None:
        min_ = real_range[0]
    if max_ is None:
        max_ = real_range[1]

    return min_, max_

def histogram(numbers, bins, range_= None, calculate_freqs=False,
              remove_outliers=False):
    '''An alternative implementation to the numpy.histogram.

    The main difference is that it accepts iterators
    '''
    min_, max_ = _calculate_range_for_histogram(numbers, range_,
                                                remove_outliers)

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
    len_numbers = 0
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
        len_numbers += 1

    if calculate_freqs:
        len_numbers = float(len_numbers)
        distrib = [(num / bin_span) /len_numbers for num in distrib]

    return (distrib, bin_edges)

def draw_histogram(values, bin_edges, title=None, xlabel= None, ylabel=None,
                   fhand=None):
    'It draws an histogram and if the fhand is given it saves it'
    modules = sys.modules.keys()
    if 'matplotlib' not in modules:
        import matplotlib
        if fhand:
            matplotlib.use('AGG')
    if 'matplotlib.pyplot' not in modules:
        from matplotlib import pyplot as plt
    else:
        plt = sys.modules['matplotlib.pyplot']

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
    modules = sys.modules.keys()
    if 'matplotlib' not in modules:
        import matplotlib
        if fhand:
            matplotlib.use('AGG')
    if 'matplotlib.pyplot' not in modules:
        from matplotlib import pyplot as plt
    else:
        plt = sys.modules['matplotlib.pyplot']

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

def _float_to_str(num, precision=2):
    'It returns a float representation'
    if 1 < num < 100:
        format_ = '%%.%if' % precision
    else:
        format_ = '%%.%ie' % precision
    return format_ % num

def _calculate_boxplot_percentiles(numpy_vects, stats_fhand, xlabels=None):
    'It calculates the boxplot stats from the given lists'

    sep = '\t'

    stats_fhand.write(sep.join(['distrib', 'mean', 'std_deviation',
                                '1st_quartile', 'median', '3rd_qualtile']))
    stats_fhand.write('\n')
    for index, vect in enumerate(numpy_vects):
        if xlabels:
            index = int(xlabels[index])
        result = ['%.2i' % index]
        result.append(_float_to_str(numpy.mean(vect)))
        result.append(_float_to_str(numpy.std(vect)))
        #the percentiles are calculated with. The vector may be empty, so we
        #don't calculate the percentile
        if vect.any():
            quarts = _calculate_percentiles(vect, [25, 50, 75])
            result.extend([_float_to_str(num) for num in quarts])

        stats_fhand.write(sep.join(result))
        stats_fhand.write('\n')

def draw_boxplot(vectors_list, fhand=None, title=None, xlabel= None,
                 ylabel=None, stats_fhand=None, xlabels=None,
                 max_plotted_boxes=None):
    'Given a list of lists it draws a boxplot'

    if not vectors_list:# or not vectors_list[0]:
        raise ValueError('No values to process')

    modules = sys.modules.keys()
    if 'matplotlib' not in modules:
        import matplotlib
        if fhand:
            matplotlib.use('AGG')
    if 'matplotlib.pyplot' not in modules:
        from matplotlib import pyplot as plt
    else:
        plt = sys.modules['matplotlib.pyplot']

    numpy_vects = [numpy.ravel(vect) for vect in vectors_list]

    if stats_fhand:
        _calculate_boxplot_percentiles(numpy_vects, stats_fhand, xlabels)

    fig = plt.figure()
    axes = fig.add_subplot(111)

    if xlabel:
        axes.set_xlabel(xlabel)
    if ylabel:
        axes.set_ylabel(ylabel)
    if title:
        axes.set_title(title)

    if max_plotted_boxes:
        step = len(numpy_vects)//max_plotted_boxes
        numpy_vects = numpy_vects[::step]
    axes.boxplot(numpy_vects)

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

def _calculate_percentiles(numbers, percents):
    'It calculates the percentiles for some numbers'
    #we need a numpy array
    if 'any' not in dir(numbers):
        numbers = numpy.ravel(numbers)
    if not numbers.any():
        raise ValueError('No data to calculate percentiles')

    modules = sys.modules
    if 'matplotlib' not in modules.keys():
        import matplotlib
        matplotlib.use('AGG')

    if 'matplotlib.mlab' not in modules.keys():
        import matplotlib.mlab as mlab
    else:
        mlab = modules['matplotlib.mlab']

    percentiles = mlab.prctile(numbers, percents)
    return list(percentiles)

class CachedArray(object):
    '''It stores numbers for future use.

    It can work in memory or write it contents to a file.
    '''
    def __init__(self, typecode):
        'The init '
        self._typecode = typecode
        self._array = array(typecode)
        self._sample = array(typecode)  #this sample will always be in memory
        self._sample_length = SAMPLE_LENGTH
        self._cache_fhand = None
        self._max = None
        self._min = None
        self._len = 0
        self._sum = 0
        self._check_in_cycles = MEMORY_CHECK_CYCLES

    def _get_sample_length(self):
        'it returns the sample length'
        return self._sample_length
    sample_length = property(_get_sample_length)

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

        #the sample
        sample = self._sample
        sample_length = self._sample_length
        if len(sample) < sample_length:
            sample.append(item)
        else:
            random_draw = random.randint(0, self._len - 2)
            if random_draw < sample_length:
                sample[random_draw] = item

    def _get_sample(self):
        'It returns the sample'
        return self._sample
    sample = property(_get_sample)

    def get_sample_percentiles(self, percents):
        'It returns the percentiles given a percent list'
        if not self._sample:
            raise ValueError('No data to calculate percentiles')

        modules = sys.modules
        if 'matplotlib' not in modules.keys():
            import matplotlib
            matplotlib.use('AGG')

        if 'matplotlib.mlab' not in modules.keys():
            import matplotlib.mlab as mlab
        else:
            mlab = modules['matplotlib.mlab']

        vect = numpy.ravel(self.sample)
        percentiles = mlab.prctile(vect, percents)
        return list(percentiles)

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

def _calculate_bins(values):
    'It calculates how many binds to use giving a list of values'
    num_values = len(values)
    if num_values < 1000:
        return ValueError('Need al least 1000 values')

    bins = int((num_values * 0.1) / 100)
    if bins < 20:
        bins = 20
    return bins

#taken from http://code.activestate.com/recipes/278260/
def _find_index(sorted_list, value, index_buffer=0):
    ''' Given a sorted_list and value, return the index i where
    sorted_list[i-1] <= value < sorted_list[i]
    Which means,
    sorted_list.insert( findIndex(sorted_list, value), value )
    will give a sorted list
    '''
    if len(sorted_list) == 2:
        if value == sorted_list[-1]:
            return index_buffer + 2
        elif value >= sorted_list[0]:
            return index_buffer + 1
        else:
            return 0
    else:
        length_ = len(sorted_list)
        first_half  = sorted_list[:length_//2 + 1]
        second_half = sorted_list[(length_//2):]

        if second_half[-1] <= value:
            return index_buffer + len(sorted_list)
        elif value < first_half[0]:
            return index_buffer
        else:
            if first_half[-1] < value:
                return _find_index(second_half, value,
                                   index_buffer=length_//2 + index_buffer)
            else:
                return _find_index(first_half, value, index_buffer=index_buffer)
