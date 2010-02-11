#!/usr/bin/env python
'''It calculates some statistics for a contig file:

    - consensus_length
    - reads_length
    - coverage_frec
'''

from optparse import OptionParser
from franklin.contig_io import get_parser
#from franklin.SeqVariation import contig_to_read_list
import matplotlib.pyplot as plt

def contig_stats(contig, data_dict):
    '''It takes the contig and calculate the stats '''
    if  data_dict is None:
        data_dict = {}
        data_dict['consensus_length'] = 0
        data_dict['reads_length']     = 0
        data_dict['coverage_frec']       = {}
    #consensus length
    consensus_len = len(contig.consensus.sequence)
    data_dict['consensus_length'] += consensus_len
    #reads total length
    for seq in contig:
        data_dict['reads_length'] += len(seq)
    #coverage_frec
    data_dict['coverage_frec'] = coverage_frec(contig,
                                               data_dict['coverage_frec'])

    return  data_dict

def longest_read(alignment):
    ''' It returns the longest string length in the list'''
    longest = 0
    for read in alignment:
        len_read = len(read)
        if len_read > longest:
            longest = len_read
    return longest

def coverage_frec(contig, coverage_data):
    '''It calculates the coverafe frecuency and it puts into the data dict '''
    #proxycontig = contig_to_read_list(contig)
    proxycontig = contig
    longest_read_len = longest_read(proxycontig)
    for i in range(longest_read_len):
        cont = 0
        for seq in proxycontig:
            try:
                nucleotide = seq[i]
            except IndexError:
                continue
            if nucleotide.isalpha():
                cont += 1
        if cont != 0:
            if cont in coverage_data:
                coverage_data[cont] += 1
            else:
                coverage_data[cont] = 1
    return coverage_data

def dict2list(coverage_data):
    '''It retusn a list with the dictionary keys sorted by its key '''
    data = []
    for i in range(max(coverage_data.keys())):
        if i in coverage_data:
            data.append(coverage_data[i])
        else:
            data.append(0)
    return data

def plot_frec(coverage_data):
    '''It plots the coverage frecuency data'''
    data_list = dict2list(coverage_data)
    plt.bar(left=range(max(coverage_data.keys())), height=data_list,
            align='center', color='y', width=0.4)
    plt.ylabel("Num of nucleotides coveraged")
    plt.title("Nucleotide coverage frecuency")
    plt.show()

def main():
    '''The main section'''
    parser = OptionParser('usage: %prog [-v] -nNODES ...', version='%prog 1.0')
    parser.add_option('-i', '--infile', dest='infile',
                      help='Caf/Ace input file')
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile

    print "Starting file indexing"
    parser = get_parser(open(infile, 'r'), infile[-3:].lower() )
    print "File indexing finished"

    data_dict = None
    for contig in parser.contigs():
        data_dict    = contig_stats(contig, data_dict)

    print data_dict['consensus_length']
    print data_dict['reads_length']
    print data_dict['coverage_frec']
    plot_frec(data_dict['coverage_frec'])

if __name__ == '__main__':
    main()
