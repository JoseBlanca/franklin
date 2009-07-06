'''This module holds utilities to filter sequences.

The filtering can be divided in two kinds of functions, the ones that return
a bool meaning if the sequence should be removed or not and the ones that
return a masked sequence or trimmed sequence. In this latter case the filter
can also return None if no sequence is left after the filtering process.
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

from itertools import imap
from StringIO  import StringIO
import re

#from biolib.biolib_utils import temp_fasta_file, NamedTemporaryDir, parse_fasta
from biolib.biolib_cmd_utils import  create_runner
from biolib.alignment_search_result import (ExonerateParser, BlastParser,
                                            FilteredAlignmentResults)

#def ifiltering_map(func, *iterators):
#    """Version of imap that only yield the True items."""
#    for result in imap(func, *iterators):
#        if result:
#            yield result
#pylint: disable-msg=R0903
#class BlastRunner(object):
#    'A blast runner for sequence objects'
#    blast_param_def = {'database':{'required':True,  'option': '-d'},
#                       'program' :{'required':True,  'option':'-p'},
#                       'expect'  :{'default': 0.0001,'option': '-e'},
#                       'nhitsv'  :{'default': 20,    'option':'-v'},
#                       'nhitsb'  :{'default': 20,    'option':'-b'},
#                       'megablast':{'default':'T',  'option':'-n'}
#                       }
#
#    def __init__(self, parameters):
#        '''The init.
#
#        Those are the parameters that are required by blast
#        database -- The blast database name
#        program  -- The blast program (e.g blastn or blastp)
#
#        These parameters are required by us, so we need to fill them
#        expect   -- The expect to filter the result
#        nhits    -- number of results to keep (default 20)
#        use_megablast -- default(False)
#        '''
#
#        why don't we just use a run_blast function?
#        because this interface would also work with if we already have a
#        multiblast result file and we want to access it in a random way.
#
#        if 'database' not in parameters or parameters['database'] == '':
#            raise RuntimeError('database parameter is required')
#        else:
#            self._database = parameters['database']
#        if 'program' not in  parameters or parameters['program'] == '':
#            raise RuntimeError('blast program is required')
#        else:
#            self._program  = parameters['program']
#        if 'expect' not in parameters or parameters['expect'] == '':
#            raise RuntimeError('expect parameter is required to create filters')
#        else:
#            self._expect  = parameters['expect']
#        if 'nhits' not in  parameters or parameters['nhits'] == '':
#            self._nhits = 20
#        else:
#            self._nhits = parameters['nhits']
#        if 'use_megablast' not in  parameters or \
#                                            parameters['use_megablast'] == '':
#            self._use_megablast = False
#        else:
#            self._use_megablast = True
#
#    def get_result(self, sequence):
#        'Given a sequence it returns the xml blast result as a string'
#        we create the fasta file
#        fastah = temp_fasta_file(sequence)
#        we run the blast
#        nhits = str(self._nhits)
#        cmd = ['blastall', '-i', fastah.name, '-p', self._program,
#               '-e', str(self._expect), '-m', '7', '-v', nhits, '-b', nhits,
#               '-d', self._database]
#        if self._use_megablast:
#            cmd.append('-n')
#            cmd.append('T')
#        stdout, stderr, retcode = call(cmd)
#        if retcode:
#            raise RuntimeError('Problem running blastall: ' + stderr)
#        fastah.close()
#        return stdout
#
#def create_blast_filter(parameters , keep_better_hits=True):
#    '''A function factory factory that creates blast filters.
#
#    It returns a function that will accept a sequence and it will return
#    True or False depending on the blast outcome.
#    database is the blast database name.
#    program is the blast program (e.g. blastn or blastp)
#    '''
#    xml_blast_source = BlastRunner(parameters)
#    def blast_filter(sequence):
#        'Given a sequence it returns True or False depending on the blast'
#        first we need the xml blast result
#        xml_blast = StringIO(xml_blast_source.get_result(sequence))
#        now we want a summary
#        biopython_blast = NCBIXML.parse(xml_blast)
#        summary = BlastSummary(biopython_blast.next())
#        we filter it
#        expect = parameters['expect']
#        summary.filter_expect_threshold(expect)
#        is the filter positive or not?
#        if len(summary.hits):
#            result = True
#        else:
#            result = False
#        if not keep_better_hits:
#            result = not result
#        return result
#    return blast_filter
#
#pylint: disable-msg=R0903
#class ExonerateRunner():
#    '''A exonerate runner for sequences'''
#    def __init__(self, parameters):#target, model=None, show_options=None):
#        ''' Initiator. Here you defines how is going to run your exonerate
#         instance
#
#         By default we are going to use an gaped alignment model. We will use
#          affine:local , similar to the Smith-Waterman-Gotoh algorithm.
#         Show options, parameter is optional as well. By default we will use
#         a format easy to arse for us
#         '''
#        if 'target' not in parameters or parameters['target'] == '':
#            raise RuntimeError('target parameter is required')
#        else:
#            self._target = ['--target', parameters['target']]
#        if 'model' not in parameters or parameters['model'] == '':
#            self._model = [ '--model', 'affine:local']
#        else:
#            self._model = ['--model', parameters['model']]
#        if 'show_option' not in parameters or parameters['show_options'] == '':
#            self._show_options = ['--showalignment', 'False', '--showvulgar',
#                                'False', '--ryo', "cigar_like:%S %ql %tl\n"]
#        else:
#            self._show_options = parameters['show_options']
#
#    def get_result(self, sequence):
#        '''Giving a sequence in returns the exonerate result with the show
#        option you have decide'''
#        fastah = temp_fasta_file(sequence)
#        query = ['--query', fastah.name]
#
#        cmd = ['exonerate']
#        cmd.extend(self._model)
#        cmd.extend(self._show_options)
#        cmd.extend(self._target)
#        cmd.extend(query)
#        stdout, stderr, retcode = call(cmd)
#        if retcode:
#            raise RuntimeError('Problem running exonerate: ' + stderr)
#        fastah.close()
#        return stdout
#SSAHA2_OPTIONS = {'adaptors':{'builder': ['-kmer', '4'],
#                              'ssaha': ['-seeds', '2', '-score', '10',
#                                        '-sense', '1', '-cmatch', '10',
#                                        '-ckmer', '6', '-identity', '90',
#                                        '-depth', '5', '-cut', '999999999',
#                                        '-memory', '500', ]}
#                 }
#class SsahaRunner(object):
#    'It creates a ssaha2 runner.'
#
#    def __init__(self, subject, options=None):
#        # We have to change this to the new model with parameters instead
#         subject/options
#        '''The init
#
#        The subject can be a string or a file with fasta sequences, do not use
#        StringIO. It  will be hashed using ssaha2Build.
#        The options should be a dict with two list of parameters one the hash
#        table builder and other for ssaha itself.
#        Other way to give the options is a string with the name of a
#        precompiled option like 'adaptors'.
#        '''
#        the options
#        if options is None:
#            options = {'builder':[], 'ssaha':[]}
#        if not isinstance(options, dict):
#            it should be a string
#            if options in SSAHA2_OPTIONS:
#                options = SSAHA2_OPTIONS[options]
#            else:
#                raise ValueError('No precompiled options for options: ' +
#                                 options)
#        self._options = options
#        self._base_hash_file = None
#        self._create_hash_file(subject)
#
#    def _create_hash_file(self, subject):
#        '''It creates a hash file using ssaha2Build.
#
#        It requires a subject with a fileh to a fasta file or a string with
#        the fasta sequences.
#        '''
#        if 'name' not in dir(subject):
#            subject is not a file, we have to create one
#            tempf = tempfile.NamedTemporaryFile(prefix='subject_seq',
#                                                suffix='.fasta')
#            tempf.write(subject)
#            tempf.flush()
#            subject = tempf
#        now we can create the hash table
#        hash_file = tempfile.NamedTemporaryFile()
#        cmd = ['ssaha2Build']
#        the options
#        cmd.extend(self._options['builder'])
#        cmd.extend(['-save', hash_file.name, subject.name])
#        print ' '.join(cmd)
#        pylint: disable-msg=W0612
#        try:
#            stdout, stderr, retcode = call(cmd)
#        except OSError:
#            raise OSError('You have not properly instaled ssaha2Build')
#        we don't need the stdout
#        if retcode:
#            raise RuntimeError('Problem running ssaha2Build: ' + stderr)
#        self._base_hash_file = hash_file
#        print "hl"
#
#    def __del__(self):
#        'We have to remove some files associated with the hash'
#        base = self._base_hash_file
#        for extension in ('head', 'body', 'name', 'base', 'size'):
#            fpath = base.name + '.' + extension
#            if os.path.exists(fpath):
#                os.remove(fpath)
#
#    def get_result(self, sequence):
#        'Given a sequence it returns the ssaha2 output as a string'
#        we create the fasta file
#        fastah = temp_fasta_file(sequence)
#        we run the ssaha2
#        cmd = ['ssaha2']
#        the options
#        cmd.extend(self._options['ssaha'])
#        the hash table
#        cmd.append('-save')
#        cmd.append(self._base_hash_file.name)
#        the query file
#        cmd.append(fastah.name)
#        try:
#            stdout, stderr, retcode = call(cmd)
#        except OSError:
#            raise OSError('You have not properly instaled ssaha2')
#        if retcode:
#            raise RuntimeError('Problem running ssaha2: ' + stderr +
#                               '\ncommand was: ' + ' '.join(cmd))
#        fastah.close()
#        return stdout

def create_filter(aligner_cmd, cmd_parameters, match_filters=None,
                            result_filters=None):
    '''A function factory factory that creates exonerate filters.

    It returns a function that will accept a sequence and it will return
    True or False depending on the exonerate outcome.
    parameters is a dictionary and key are defined in ExonerateRunner.
    Required is only the target fasta file
    '''
    #runners = {'blast':BlastRunner, 'exonerate':ExonerateRunner}

    parsers = {'blast':BlastParser, 'exonerate':ExonerateParser}
    parser  = parsers[aligner_cmd]
    binary  = {'blast':'blast2', 'exonerate':'exonerate'}

    Runner = create_runner(kind=aligner_cmd, bin_=binary[aligner_cmd])
    source = Runner(parameters=cmd_parameters)

    def _filter(sequence):
        'Giving a sequence it returns true or False depending on the exonerate'

        source_result    = StringIO(source.get_result(sequence))
        results          = parser(source_result)
        filtered_results = FilteredAlignmentResults(filters=match_filters,
                                                   results=results)
        try:
            result = filtered_results.next()
            if not len(result['matches']):
                return False
        except StopIteration:
            return False

        filter_result = _filtered_match_results(filters=result_filters,
                                                result=result)
        return filter_result

    return _filter


def _filtered_match_results(filters, result):
    '''It returns True or False depending if the result pass the filter or
    not'''
    if filters is None:
        return True
    num_matches = len(result['matches'])
    for filter_ in filters:
        if filter_['kind'] == 'num_matches':
            min_num_matches = filter_['value']
    if num_matches < min_num_matches:
        return False
    else:
        return True



#### Functions that strip the sequences
def strip_seq_by_quality(seq, quality_treshold, min_quality_bases=None,
                       min_seq_length=None, quality_window_width=None):
    '''It trims from the seq bad quality extremes

    '''
    quality = seq.qual
    seq_len = len(quality)
    # Define defaults
    (min_quality_bases, min_seq_length, quality_window_width) = \
    _trim_seq_by_quality_defaults(seq_len, min_quality_bases, min_seq_length,
                                  quality_window_width)
    #calculate the sliding window qual
    qual_average_array = _calculate_sliding_window_qual(quality,
                                                        quality_window_width)
    # Now we convert the quality average array in a boolean string. 1 if the
    # quality is good, and 0 if it is bad. Depending on the treshold
    boolean_quality_treshold = _quality_to_boolean(quality_treshold,
                                                   qual_average_array)

    start, end = _trim_bad_qual_extremes(boolean_quality_treshold,
                                         min_quality_bases)
    new_seq = seq[start:end]
    if len(new_seq.qual) < min_seq_length:
        return None
    else:
        return new_seq
def _trim_bad_qual_extremes(bool_seq, min_quality_bases):
    '''It returns start and and of the new '''
    start = re.search('0*(1{%d})' % min_quality_bases, bool_seq).start(1)

    #First I reverse de quality string
    bool_seq = bool_seq[::-1]
    end = re.match('0*(1{%d})' % min_quality_bases , bool_seq).start(1)
    end = len(bool_seq) - end
    return start, end
def _calculate_sliding_window_qual(quality, quality_window_width):
    '''It takes into account the neighbour qualities to calculate the new
    quality. The window size is variable '''
    average_quality = []
    #pylint: disable-msg=W0612
    for i, single_qual in enumerate(quality):
        if quality_window_width > i:
            qual_slice = quality[0: i + quality_window_width +1]
        else:
            qual_slice = quality[i -quality_window_width:
                                 i + quality_window_width +1]
        average = sum(qual_slice) / float(len(qual_slice))
        average_quality.append(average)
    return average_quality
def _quality_to_boolean(quality_treshold, qual_average_array):
    '''It converts the average_quality array in a boolean string.  1 if the
    quality is good, and 0 if it is bad. Depending on the treshold'''
    boolean_qual = []
    for quality in qual_average_array:
        if quality > quality_treshold:
            boolean_qual.append('1')
        else:
            boolean_qual.append('0')
    return ''.join(boolean_qual)
def _trim_seq_by_quality_defaults(seq_len, min_quality_bases, min_seq_length,
                                  quality_window_width):
    '''We define defaults for the  trim_seq_by_quality function'''
    if min_quality_bases is None:
        if seq_len < 70:
            min_quality_bases = 6
        elif seq_len < 300:
            min_quality_bases = 15
        else:
            min_quality_bases = 2

    if quality_window_width is None:
        if seq_len < 70:
            quality_window_width = 2
        elif seq_len < 300:
            quality_window_width = 3
        else:
            quality_window_width = 4

    if min_seq_length is None:
        fiftypercent = seq_len / 2
        if fiftypercent > 100:
            min_seq_length = 100
        else:
            min_seq_length = fiftypercent
    return (min_quality_bases, min_seq_length, quality_window_width)


