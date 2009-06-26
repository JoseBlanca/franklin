'''
Created on 2009 mai 13

@author: peio
'''
from optparse import OptionParser
from biolib.contig_io import get_parser_by_name
from biolib.biolib_utils import fasta_str
import sys

def main():
    '''The main function '''
    parser = OptionParser('usage: %prog -i [caffile|acefile] -o fastafile')
    parser.add_option('-o', '--outfile', dest='outfile', help='Output file')
    parser.add_option('-i', '--infile', dest='infile', help='Maximun limit')
    parser.add_option('-c', '--contigs', dest='contig_list', \
                      help='Contig list to print. Comma separated')
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile       = options.infile
    if options.outfile is None:
        outfile = sys.stdout
    else:
        outfile = options.outfile
        outfile = open(outfile, 'w')

    if options.contig_list is None:
        contig_list = None
    else:
        contig_list = options.contig_list.split(',')

    parser = get_parser_by_name(infile)
    for contig in parser.contigs():
        name = contig.consensus.sequence.name
        if contig_list is  None or name in contig_list:
            sequence = contig.consensus.sequence
            outfile.write(fasta_str(sequence, name))

    outfile.close()

if __name__ == '__main__':
    main()
