'''
Created on 2009 api 7

@author: peio
'''
from biolib.SeqVariation import seqvariations_in_alignment
from biolib.contig_parser import CafParser, AceParser
#from sys import argv 

def main():
    ''' Main function where we find snps'''
    #caf_file = argv[1]
    #caf_file = "/home/peio/work_in/assemblers/eucalyptus_c5476.caf"
    #caf_file = "/home/peio/work_in/assemblers/eucalyptus_out.caf"
    #caf_file = '/home/jose/tmp/example.caf'
    caf_file = '/home/peio/devel/biolib/src/data/example3.caf'
    ace_file = '/home/peio/devel/biolib/src/data/example3.ace'
    print "starting file indexing"
    parser = AceParser(ace_file)
    parser = CafParser(caf_file)
    
    print "file indexing finished"
    
    var_count = {}
    contig_count = 0
    for contig in parser.contigs():
        contig_count += 1
        for seqvar in seqvariations_in_alignment(contig):
            kind = seqvar.kind()
            print seqvar.location
            if not kind in var_count:
                var_count[kind] = 0
            var_count[kind] += 1
        print 'contigs: ' + str(contig_count)
        print 'variations:' + str(var_count)
        break


if __name__ == '__main__':
    import cProfile
    #cProfile.run('main()')
    main()

