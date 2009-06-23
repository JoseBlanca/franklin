'''This module holds utilities to filter sequences.

The filtering can be divided in two kinds of functions, the ones that return
a bool meaning if the sequence should be removed or not and the ones that
return a masked sequence or trimmed sequence. In this latter case the filter
can also return None if no sequence is left after the filtering process.
'''

from itertools import imap
from StringIO  import StringIO
import os, tempfile

from biolib.biolib_utils import call, temp_fasta_file
from biolib.blast_summary import BlastSummary
from Bio.Blast import NCBIXML

def ifiltering_map(func, *iterators):
    """Version of imap that only yield the True items."""
    for result in imap(func, *iterators):
        if result:
            yield result

class BlastRunner(object):
    'A blast runner for sequence objects'
    def __init__(self, database, program, expect, nhits=20,
                 use_megablast=False):
        '''The init.
        
        keyword arguments
        database -- The blast database name
        program  -- The blast program (e.g blastn or blastp)
        expect   -- The expect to filter the result
        nhits    -- number of results to keep (default 20)
        use_megablast -- default(False)
        '''
        #pylint: disable-msg=R0913
        #blast requires this arguments, and maybe some more

        #why don't we just use a run_blast function?
        #because this interface would also work with if we already have a 
        #multiblast result file and we want to access it in a random way.
        self._database = database
        self._program  = program
        self._expect   = expect
        self._nhits    = nhits
        self._use_megablast = use_megablast

    def get_blast(self, sequence):
        'Given a sequence it returns the xml blast result as a string'
        #we create the fasta file
        fastah = temp_fasta_file(sequence)
        #we run the blast
        nhits = str(self._nhits)
        cmd = ['blastall', '-i', fastah.name, '-p', self._program,
               '-e', str(self._expect), '-m', '7', '-v', nhits, '-b', nhits,
               '-d', self._database]
        if self._use_megablast:
            cmd.append('-n')
            cmd.append('T')
        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('Problem running blastall: ' + stderr)
        fastah.close()
        return stdout

def create_blast_filter(expect, database, program, keep_better_hits=True,
                        use_megablast=True):
    '''A function factory factory that creates blast filters.
    
    It returns a function that will accept a sequence and it will return
    True or False depending on the blast outcome.
    database is the blast database name.
    program is the blast program (e.g. blastn or blastp)
    '''
    xml_blast_source = BlastRunner(database=database, program=program,
                                    expect=expect, use_megablast=use_megablast)
    def blast_filter(sequence):
        'Given a sequence it returns True or False depending on the blast'
        #first we need the xml blast result
        xml_blast = StringIO(xml_blast_source.get_blast(sequence))
        #now we want a summary
        biopython_blast = NCBIXML.parse(xml_blast)
        summary = BlastSummary(biopython_blast.next())
        #we filter it
        summary.filter_expect_threshold(expect)
        #is the filter positive or not?
        if len(summary.hits):
            result = True
        else:
            result = False
        if not keep_better_hits:
            result = not result
        return result
    return blast_filter

SSAHA2_OPTIONS = {'adaptors':{'builder': ['-kmer', '6'],
                              'ssaha': ['-seeds', '2', '-score', '15',
                                        '-sense', '1', '-cmatch', '10',
                                        '-ckmer', '6', '-identity', '90',
                                        '-depth', '5', '-cut', '999999999',
                                        '-memory', '500']}
                 }

class SsahaRunner(object):
    'It creates a ssaha2 runner.'

    def __init__(self, subject, options=None):
        '''The init

        The subject can be a string or a file with fasta sequences, do not use
        StringIO. It  will be hashed using ssaha2Build.
        The options should be a dict with two list of parameters one the hash
        table builder and other for ssaha itself.
        Other way to give the options is a string with the name of a
        precompiled option like 'adaptors'.
        '''
        #the options
        if options is None:
            options = {'builder':[], 'ssaha':[]}
        if not isinstance(options, dict):
            #it should be a string
            if options in SSAHA2_OPTIONS:
                options = SSAHA2_OPTIONS[options]
            else:
                raise ValueError('No precompiled options for options: ' +
                                 options)
        self._options = options
        self._base_hash_file = None
        self._create_hash_file(subject)

    def _create_hash_file(self, subject):
        '''It creates a hash file using ssaha2Build.
        
        It requires a subject with a fileh to a fasta file or a string with
        the fasta sequences.
        '''
        if 'name' not in dir(subject):
            #subject is not a file, we have to create one
            tempf = tempfile.NamedTemporaryFile(prefix='subject_seq',
                                                suffix='.fasta')
            tempf.write(subject)
            tempf.flush()
            subject = tempf
        #now we can create the hash table
        hash_file = tempfile.NamedTemporaryFile()
        cmd = ['ssaha2Build']
        #the options
        cmd.extend(self._options['builder'])
        cmd.extend(['-save', hash_file.name, subject.name])
        print ' '.join(cmd)
        #pylint: disable-msg=W0612
        stdout, stderr, retcode = call(cmd)
        #we don't need the stdout
        if retcode:
            raise RuntimeError('Problem running ssaha2Build: ' + stderr)
        self._base_hash_file = hash_file

    def __del__(self):
        'We have to remove some files associated with the hash'
        base = self._base_hash_file
        for extension in ('head', 'body', 'name', 'base', 'size'):
            fpath = base.name + '.' + extension
            if os.path.exists(fpath):
                os.remove(fpath)

    def get(self, sequence):
        'Given a sequence it returns the ssaha2 output as a string'
        #we create the fasta file
        fastah = temp_fasta_file(sequence)
        #we run the ssaha2
        cmd = ['ssaha2']
        #the options
        cmd.extend(self._options['ssaha'])
        #the hash table
        cmd.append('-save')
        cmd.append(self._base_hash_file.name)
        #the query file
        cmd.append(fastah.name)
        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('Problem running ssaha2: ' + stderr + 
                               '\ncommand was: ' + ' '.join(cmd))
        fastah.close()
        return stdout

def create_ssaha_filter(subject, options=None):
    '''This function factory returns a ssaha filter for sequences.

    The subject can be a string or a file with fasta sequences, do not use
    StringIO. It  will be hashed using ssaha2Build.
    The options should be a list with parameters for ssaha or a string
    with the name of a precompiled option like 'adaptor'.
    '''
    ssaha_source = SsahaRunner(subject=subject, options=options)

    def ssaha_filter(sequence):
        'Given a sequence it returns True or False depending on the ssaha'
        #first we need the ssaha result
        ssaha = ssaha_source.get(sequence)
        print ssaha
    return ssaha_filter
