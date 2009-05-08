#!/usr/bin/env python

from optparse import OptionParser
from sys import stdin
import matplotlib.pyplot as plt

def main():
    ''' Main function where we find snps'''
    
    parser = OptionParser('usage: %prog [-v] -nNODES ...', version='%prog 1.0')
    parser.add_option('-i', '--infile', dest='infile',
                      help='Input file')
    parser.add_option('-o', '--outfile',  dest='outfile', 
                      help='Output file')
    parser.add_option('-M', '--max',  dest='max',  type="float",
                      help='Maximun limit')
    parser.add_option('-m', '--min', dest='min',type="float",
                      help='Minimun Limit')
    parser.add_option('-t', '--interval', dest = 'interval', type='int',  
                      help = 'plot interval?'  )
    (options, args) = parser.parse_args()
    
    
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
    n, bins, patches = axes.hist(data_list, bins=interval, range=plot_range, 
                                 facecolor='green', alpha=0.75)
#    axes.set_xlabel('Second allele ')
#    axes.set_ylabel('Num of ocurrences')
    plt.show()

    
    
if __name__ == '__main__':
    main()