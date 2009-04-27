'''
Created on 2009 api 7

@author: peio
'''
from biolib.SeqVariation import seqvariations_in_alignment, contig_to_read_list
from biolib.contig_parser import CafParser, AceParser
from biolib.contig_cleaner import contig_strip
#from sys import argv 

def main():
    ''' Main function where we find snps'''
    #caf_file = argv[1]
    #caf_file = "/home/peio/work_in/assemblers/eucalyptus_c5476.caf"
#    caf_file = "/home/peio/work_in/assemblers/eucalyptus_out.caf"
    #caf_file = '/home/jose/tmp/example.caf'
    caf_file = '/home/peio/devel/biolib/src/data/example3.caf'
    ace_file = '/home/peio/devel/biolib/src/data/example3.ace'
    print "starting file indexing"
    parser = AceParser(ace_file)
    parser = CafParser(caf_file)
    
    print "file indexing finished"
    
    contig_count = 0
    for contig in parser.contigs():
        var_count = {}
        contig_count += 1
#        for read in contig:
#            print read.location
#            print read.location.forward
        print contig
        contig = contig_strip(contig, 5)
#        print contig
            
        print "Searching snps"
        for seqvar in seqvariations_in_alignment(contig):
            kind = seqvar.kind()
            if not kind in var_count:
                var_count[kind] = 0
            var_count[kind] += 1
        
        print 'variations:', var_count
    print 'contigs: ' , contig_count

if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()

