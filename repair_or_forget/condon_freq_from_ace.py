#!/usr/bin/env python
'''
It returns a matrix with the codon frecuency, giving a caf file with est
'''
from optparse import OptionParser
from biolib.contig_io import get_parser
from biolib.biolib_cmd_utils import get_best_orf

def main():
    '''The main fucntion '''
    parser = OptionParser('usage: %prog [-v] -nNODES ...', version='%prog 1.0')
    parser.add_option('-M', '--matrix', dest='matrix',
                      help='codon frecuency matrix')
    parser.add_option('-o', '--outfile', dest='outfile', help='Output file')
    parser.add_option('-i', '--infile', dest='infile', help='Maximun limit')
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile
    if options.matrix is None:
        parser.error('ESTScan need a codon frecuency matrix')
    else:
        matrix = options.matrix

    print "Starting file indexing"
    parser = get_parser(infile, infile[-3:].lower())
    print "File indexing finished"

    codons = {}
    for contig in parser.contigs():
        sequence = contig.consensus.sequence
        orf_dna = get_best_orf(sequence, matrix_path=matrix)[0]
        cont = 0

        while True:
            codon = orf_dna[cont:cont + 3]
            codon = codon.upper()
            if len(codon) != 3:
                break
            if 'X' not in codon:
                if codon not in codons:
                    codons[codon] = 1
                else:
                    codons[codon] += 1
            cont += 3

    print codons

if __name__ == '__main__':
    main()
