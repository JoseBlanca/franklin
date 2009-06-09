# -*- coding= UTF-8
'''
Created on 2009 eka 2

@author: peio
'''
from optparse import OptionParser
from biolib.contig_parser import CafParser, AceParser
from biolib.SeqVariation import contig_to_read_list, longest_read
import matplotlib.pyplot as plt

def contig_stats(contig, data_dict):
    '''It takes the contig and calculate teh stats '''
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

def coverage_frec(contig, coverage_data):
    '''It calculates the coverafe frecuency and it puts into the data dict '''
    proxycontig = contig_to_read_list(contig)
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
    parser.add_option('-i', '--infile', dest='infile', help='Maximun limit')
    (options, args) = parser.parse_args()

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile
    
    print "Starting file indexing"  
    if infile[-3:].lower() == 'ace':
        print "File type: ace"
        parser = AceParser(infile)
    elif infile[-3:].lower() == 'caf':
        print "File type: caf"
        parser = CafParser(infile)
    print "File indexing finished"
    
    data_dict = None
    for contig in parser.contigs():
        data_dict    = contig_stats(contig, data_dict)
    
    
    print data_dict['consensus_length']
    print data_dict['reads_length']
    print data_dict['coverage_frec']
    plot_frec(data_dict['coverage_frec'])
    
    # Si ya tienes el print de  tu diccionario calculado, commenta todo lo 
    # demás y pon el dict aquí para visualizarlo
    #data_dict=
    #plot_frec(data_dict)
    
if __name__ == '__main__':
    main()