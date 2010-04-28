#!/usr/bin/env python
'Given a column with numbers in a file or stdin it plots an histogram'

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

from optparse import OptionParser
from sys import stdin
import matplotlib.pyplot as plt

def main():
    ''' Main'''

    parser = OptionParser('usage: %prog [-v] -nNODES ...', version='%prog 1.0')
    parser.add_option('-i', '--infile', dest='infile',
                      help='Input file')
    parser.add_option('-o', '--outfile',  dest='outfile',
                      help='Output file')
    parser.add_option('-M', '--max',  dest='max',  type="float",
                      help='Maximun limit')
    parser.add_option('-m', '--min', dest='min', type="float",
                      help='Minimun Limit')
    parser.add_option('-t', '--interval', dest = 'interval', type='int',
                      help = 'plot interval?'  )
    options = parser.parse_args()[0]


    if options.infile is None:
        ifhand = stdin
    else:
        ifname = options.infile
        ifhand = open(ifname,'rt')
    file_content = ifhand.read()

    if options.interval is not None:
        interval = options.interval
    else:
        interval = 20

    data_list = []
    for line in file_content.split('\n'):
        line.strip()
        if line.isspace() or line == '':
            continue
        number = float(line)
        data_list.append(number)

    if options.max is None:
        opt_max = max(data_list)
    else:
        opt_max  = options.max
    if options.min is None:
        opt_min = min(data_list)
    else:
        opt_min  = options.min
    plot_range = opt_min, opt_max

    #ploting the figure
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.hist(data_list, bins=interval, range=plot_range,
                                 facecolor='green', alpha=0.75)
    plt.show()

if __name__ == '__main__':
    main()
