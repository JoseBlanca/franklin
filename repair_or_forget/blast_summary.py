'''This class '''
from Bio.Blast import NCBIXML
import re

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
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

print '''This module is DEPRECATED, port the summary_to_gff3 to the new blast summary infrastructure'''

class BlastSummary(object):
    '''It holds the most relevant information about a blast, query and subj ids,
    expects, similarities, etc.'''
    
    def __init__(self, blast):
        '''Constructor
        blast -- a biopython blast_result
        '''
        self._blast = blast
        # The summary is a list of sequences producing significant alignments
        # (hits)
        self._hits             = [] 
        self._query_length     = None
        self._query_name       = None
        self._query_definition = None
        self._summarize()
        self._blast = None # Let's free some memory

    def get_hits(self):
        '''It returns the hits for this blast'''
        return self._hits
    hits = property(get_hits)
    def get_query_length(self):
        '''It returns the query length in number of residues'''
        return self._query_length
    query_length = property(get_query_length)
    def get_query_name(self):
        '''It returns the query_id for this blast summary'''
        return self._query_name
    query_name = property(get_query_name)
    def get_query_definition(self):
        '''It returns the query definition for this blast summary'''
        return self._query_definition
    query_definition = property(get_query_definition)
    def _summarize(self):
        '''It creates the summary. The summary has the structure:
        
        hits -- An structure with the information about the query-subj pairs.
        hits is a list, each item corresponds to a sequence producing
        significant alignments. It should be sorted following the original
        blast order.
        Each item in hits is a dict with the keys:
            query_id, sub_id, hsps, query_length, subj_length.
        hsps is a list with one item for every hsp. Every item is a dict
        with the keys:
            expect (required in this case), (the following keys are not
            required for this criteria) subj_start, subj_end, query_start,
            query_end, similarity
        '''
        #the query could have a definition
        try:
            query_id, definition = self._blast.query.split(' ', 1)
        except ValueError:
            query_id   = self._blast.query
            definition = None
        
        #length of query sequence
        self._query_length     = self._blast.query_letters 
        self._query_name       = query_id
        self._query_definition = definition
        
        for alignment in self._blast.alignments:
            hit_name = alignment.title.split()[0]
            #sometimes the blast leaves an lcl| in the id
            hit_length = alignment.length
            hsps = []
            hit_start, hit_end = None, None 
            for hsp in alignment.hsps:
                expect         = hsp.expect
                subject_start  = hsp.sbjct_start
                subject_end    = hsp.sbjct_end
                query_start    = hsp.query_start
                query_end      = hsp.query_end
                hsp_length     = len(hsp.query)
                # We have to check the strand of the seq subject
                if subject_start < subject_end:
                    subject_strand = '+'
                else:
                    subject_strand = '-'
                    subject_start, subject_end = subject_end, subject_start
                # Alse the query strand
                if query_start < query_end:
                    query_strand = '+'
                else:
                    query_strand = '-'
                    query_start, query_end = query_end, query_start
                
                
                try:
                    similarity = hsp.positives*100.0/float(hsp_length)
                except TypeError:
                    similarity = None
                hsps.append({
                    'expect'         : expect,
                    'subject_start'  : subject_start,
                    'subject_end'    : subject_end,
                    'subject_strand' : subject_strand, 
                    'query_start'    : query_start,
                    'query_end'      : query_end,
                    'query_strand'   : query_strand, 
                    'similarity'   : similarity })
                # It takes the first loc and the last loc of the hsp to
                # determine hit start and end
                if hit_start is None or query_start < hit_start:
                    hit_start = query_start
                if hit_end is None or query_end > hit_end:
                    hit_end = query_end
                
            self._hits.append({
                'name'   : hit_name,
                'length' : hit_length,
                'start'  : hit_start,   
                'end'    : hit_end, 
                'expect' : hsps[0]['expect'],
                'hsps'   : hsps })
    def filter_best_expects(self, min_expect, expect_tolerance):
        '''It will remove all the query-subj pairs that doesn't meed the
        criteria. It select only the ones with the best expect values. 
        This particular filter is thought to be useful in the orthologs
        searches.
        '''
        filtered_hits = []
        from math import log10
        #the best expect correspond to the first hsps in the first sequence
        best_expect = self._hits[0]['hsps'][0]['expect']
        if best_expect == 0.0:
            log_best_expect = 0.0
        else:
            log_best_expect = log10(best_expect)
        log_tolerance = log10(expect_tolerance)
        for hit in self._hits:
            expect = hit['hsps'][0]['expect']
            if expect == 0.0:
                filtered_hits.append(hit)
                continue
            if expect >= min_expect:
                continue
            if abs(log10(expect) - log_best_expect) < log_tolerance:
                filtered_hits.append(hit)
        self._hits = filtered_hits
    def filter_expect_threshold(self, min_expect):
        '''It will remove all the query-subj pairs that doesn't meet the
        criteria. It select only the ones with the expect values above the
        threshold. 
        '''
        filtered_hits = []
        #the best expect correspond to the first hsps in the first sequence
        for hit in self._hits:
            expect = hit['hsps'][0]['expect']
            if expect <= min_expect:
                filtered_hits.append(hit)
        self._hits = filtered_hits
    def _secure_hsps(self):
        '''If there is no secure copy of the original ones it will create one'''
        import copy
        for hit in self._hits:
            if not hit.has_key('hsps_mod'):
                hit['hsps_mod'] = copy.deepcopy(hit['hsps'])
    #pylint: disable-msg=C0103
    def _filter_hits_under_similarity_threshold(self, min_similarity):
        '''It removes the hits that have no hsps above the similarity threshold.
        It also leaves a hsp_mod with the hsps above the similarity thresdhold
        '''
        self._secure_hsps() #we don't want to modify the original hsps
        filtered_hits = []
        for hit in self._hits:
            filtered_hsps = []
            for hsp in hit['hsps_mod']:
                if hsp['similarity'] >= min_similarity:
                    filtered_hsps.append(hsp)
            if len(filtered_hsps) > 0:
                hit['hsps_mod'] = filtered_hsps
                filtered_hits.append(hit)
        self._hits = filtered_hits
    #pylint: disable-msg=C0103    
    def _filter_hits_under_length_threshold(self, min_length):
        '''It will remove the hits that doesn't meet the length threshold
        '''
        self._secure_hsps() #we should use the hsp_mod copy
        filtered_hits = []
        for hit in self._hits:
            hsp_len = self._hsp_length(hit['hsps_mod'])
            if hsp_len >= min_length:
                filtered_hits.append(hit)
        self._hits = filtered_hits
    def _merge_overlaping_hsps(self):
        '''hsp 1  -------        ----->    -----------
           hsp 2       ------
           The similarity information will be lost
           The modified hsps will be in hsp_mod not in hsp
        '''
        #pylint: disable-msg=C0103
        #pylint: disable-msg=C0111
        def cmp_subj_location(x, y):
            return x['subj'] - y['subj']
        
        self._secure_hsps() #we don't want to alter the original hsps

        for hit in self._hits:
            hsp_limits = [] #all hsp starts and ends
            #we collect all start and ends
            for hsp in hit['hsps_mod']:
                hsp_limit_1 = {
                    'type' : 'start',
                    'subj' : hsp['subject_start'],
                    'query' : hsp['query_start']
                }
                hsp_limit_2 = {
                    'type' : 'end',
                    'subj' : hsp['subject_end'],
                    'query' : hsp['query_end']
                }
                hsp_limits.append(hsp_limit_1)
                hsp_limits.append(hsp_limit_2)
            #now we sort the hsp limits according their subj location
            hsp_limits.sort(cmp_subj_location)
            #now we creaet the merged hsps
            starts = 0
            merged_hsps = []
            for hsp_limit in hsp_limits:
                if hsp_limit['type'] == 'start':
                    starts += 1
                    if starts == 1:
                        subj_start = hsp_limit['subj']
                        query_start = hsp_limit['query']
                elif hsp_limit['type'] == 'end':
                    starts -= 1
                    if starts == 0:
                        subj_end = hsp_limit['subj']
                        query_end = hsp_limit['query']
                        merged_hsps.append({
                            'expect' : None,
                            'subject_start' : subj_start,
                            'subject_end' : subj_end,
                            'query_start' : query_start,
                            'query_end' : query_end,
                            'similarity' : None
                            }
                        )
            hit['hsps_mod'] = merged_hsps
    @staticmethod
    def _hsp_length(hsps):
        '''It returns the length covered by a list of hsp'''
        lenght = 0
        for hsp in hsps:
            lenght = lenght + hsp['subject_end'] - hsp['subject_start'] + 1
        return lenght
            
    #pylint: disable-msg=C0103
    def filter_similarity_threshold(self, min_length, min_similarity):
        '''It will remove all the query-subj pairs that doesn't meet the
        criteria.
        It select only the ones that share more than min_length residues
        with more than min_similarity similarity percentage. This length could
        be divided in several hsps.
        This criteria could be used to join DNA probes that could cross-
        hybridize.
        '''
        self._filter_hits_under_similarity_threshold(min_similarity)
        self._merge_overlaping_hsps()
        self._filter_hits_under_length_threshold(min_length)
    def _incompatible_length(self, hit):
        '''It returns the incompatible length. The length of the region
        that should be alignmed but it is not.
        '''
        first_hsp = hit['hsps_mod'][0]
        #subject orientation
        subj_orientation = first_hsp['subject_end'] - first_hsp['subject_start']
        if subj_orientation < 1 :
            #if sori is -1 we have a programming error. subject always is foward
            raise ValueError('Orientation error, the query should be \
                forward, but it not, this is a programming bug that should be \
                fixed')
            
        #query orientation
        query_orientation = first_hsp['query_end'] - first_hsp['query_start']
        if query_orientation > 1 :
            query_orientation = 1
        else:
            query_orientation = -1
        query_len   = self.query_length
        subject_len = hit['length']

        #First we add the first incompatible track
        #------------------------------------
        #           ||||||          ||||||
        #    ------------------------------------------
        #    *******
        if query_orientation == 1:
            #ilen -> incompatible length
            incompatible_len = first_hsp['query_start'] - 1
        else:
            incompatible_len = query_len - first_hsp['query_start']

        last_hsp = hit['hsps_mod'][-1]
        #Now we add the last incompatible track
        #                                 --- subj_last
        #------------------------------------
        #           ||||||          ||||||
        #    ------------------------------------------
        #                                 ------------- query_last
        #                                 *** to add
        subj_last = subject_len - last_hsp['subject_end']
        if query_orientation == 1:
            query_last = query_len - last_hsp['query_end']
        else:
            query_last = last_hsp['query_end'] - 1
        if subj_last <= query_last:
            to_add = subj_last
        else:
            to_add = query_last
        incompatible_len += to_add
        #Now we add the incompatible between hsps
        #------------------------------------
        #           ||||||          ||||||
        #    ------------------------------------------
        #                 **********  to add
        if len(hit['hsps']) == 1:
            return incompatible_len  #if there is only one hsp we're done
        to_add = 0

        #pylint: disable-msg=W0612
        for i, hsp in enumerate(hit['hsps_mod'][1:]):
            to_add = to_add + hit['hsps_mod'][i+1]['subject_start'] - \
                                            hit['hsps'][i]['subject_end']
        incompatible_len += to_add
        return incompatible_len
    #pylint: disable-msg=C0103   
    def _filter_hits_above_incompatibility_threshold(self, max_incompatibility):
        ''' Filter hits above imcompatible criteria'''
        self._secure_hsps()
        filtered_hits = []
        for hit in self._hits:
            #filtered_hsps = []
            incompatible_len = self._incompatible_length(hit)
            if incompatible_len < max_incompatibility:
                filtered_hits.append(hit)
        self._hits = filtered_hits
    def filter_compatibility_threshold(self, min_compatibility,
                                max_incompatibility, min_similarity):
        '''It will remove all the query-subj pairs that doesn't meet the
        criteria.
        It select only the ones that are compatible in more than
        min_compatibility. Compatible is the fragment that is aligned in an
        hsp. Also it shouldn't be more than max_incompatibility.
        Incompatibility means that they should be alignmed in that region, but
        they aren't. For instance:
            -----------------------------------
                                |||||||||||||||
                -------------------------------
                incompatible    compatible
                  region          region
        Both compatibility and incompatibility will be applied to both query
        and subject and are number of residues.
        '''
        self._filter_hits_under_similarity_threshold(min_similarity)
        self._merge_overlaping_hsps()
        self._filter_hits_under_length_threshold(min_compatibility)
        self._filter_hits_above_incompatibility_threshold(max_incompatibility)

class BlastSummaries(object):
    '''This class will yield blastSummaries when given a blast result.'''
    def __init__(self, blast):
        '''Constructor

        blast -- an XML blast file name
        '''
        self._blast_name  = blast
        blast.seek(0, 0)

        self._filters     = []
        self._blast_parse = NCBIXML.parse(blast)
        self.iter         = iter(self._blast_parse)
    def __iter__(self):
        '''It is necesary to convert our object in iterable '''
        return self
    def next(self):
        '''It returns the query and the subs for a blastSummary each time.
        
        It will return all blast results.
        It returns a tuple with the query_id and a list of subj_ids.
        If a list of filters is passed the results will be filtered following
        the given criteria.
        '''
        filters = self._filters
        summary = BlastSummary(self._blast_parse.next())
        #filters
        for filter_ in filters:
            if filter_['criteria'] == 'best_expects':
                summary.filter_best_expects(filter_['min_expect'],
                                            filter_['expect_tolerance'])
            if filter_['criteria'] == 'expect_threshold':
                summary.filter_expect_threshold(filter_['min_expect'])
            if filter_['criteria'] == 'similarity_threshold':
                summary.filter_similarity_threshold(filter_['min_length'], 
                                                    filter_['min_similarity'])
            if filter_['criteria'] == 'compatibility_threshold':
                summary.filter_compatibility_threshold(
                                                filter_['min_compatibility'],
                                                filter_['max_incompatibility'],
                                                filter_['min_similarity'])
        return summary

    def add_filter_best_expects(self, min_expect, expect_tolerance):
        '''It will remove all the query-subj pairs that doesn't meet the
        criteria. It selects only the ones with the best expect values. 
        This particular filter is thought to be useful in the orthologs
        searches.
        '''
        self._filters.append({'criteria': 'best_expects',
                              'min_expect':min_expect,
                              'expect_tolerance':expect_tolerance})
    def add_filter_expect_threshold(self, min_expect):
        '''It will remove all the query-subj pairs that doesn't meed the
        criteria. It select only the ones with the expect values above the
        threshold. 
        '''
        self._filters.append({'criteria': 'expect_threshold',
                               'min_expect':min_expect})
    def add_filter_similarity_threshold(self, min_length, min_similarity):
        '''It will remove all the query-subj pairs that doesn't meed the
        criteria.
        It select only the ones with that share more than min_length residues
        with more than min_similarity similarity percentage. This length could
        be divided in several hsps
        '''
        self._filters.append({'criteria': 'similarity_threshold',
                              'min_length':min_length,
                              'min_similarity':min_similarity})

    #pylint: disable-msg=C0103
    def add_filter_compatibility_threshold(self, min_compatibility, 
                                       max_incompatibility, min_similarity):
        '''It will remove all the query-subj pairs that doesn't meet the
        criteria.
        It select only the ones that are compatible in more than
        min_compatibility. Compatible is the fragment that is aligned in an
        hsp. Also it shouldn't be more than max_incompatibility.
        Incompatibility means that they should be alignmed in that region, but
        they aren't. For instance:
            -----------------------------------
                                |||||||||||||||
                -------------------------------
                incompatible    compatible
                  region          region
        Both compatibility and incompatibility will be applied to both query
        and subject and are percentages.
        '''
        self._filters.append({'criteria': 'compatibility_threshold',
                'min_compatibility':min_compatibility,
                'max_incompatibility':max_incompatibility,
                'min_similarity':min_similarity})

def _transform_name_with_regexp(name, regexp):
    '''It returns a matching patter of the name if it matches the regexp '''
    if regexp is None:
        return name
    else:
        reg    = re.compile(regexp)
        try:
            result = reg.match(name)
            return result.groups()[0]
        except AttributeError:
            print("regexp doesn't match")
            return name

def _get_dbxref_name(database, name):
    '''It takes db and name and returns dbxref section of the 9th column 
    of the gff3 file '''
    if database:
        return ",Dbxref=%s" % _add_db_tag(name, database)
    else:
        return  ""
def _add_db_tag(name, database):
    '''It adds the database tag to the name '''
    if database:
        return "%s:%s" % (database, name)
    else:
        return name
def summarie_to_gff3(summary, feature_type, query_db=None, subject_db=None,
                      query_regex=None, subject_regex=None):
    ''' It convert a  blast record in a gff3 output. There are arguments that 
    need some explication:
         .- feature_type : This is mandatory. And is the SOFA feature that 
            describes better the query secuence. This is because we have to put
            a feature where the other match and match_part "hang on"
        .- *_regexp: These regular expresion are use to transform the name into
            our naming schema. his is howto use it. You need to put parentesis 
            in the matching region that you need to conserve in the name.
            Example.SGD-12345 To conserve the number part you need => "SGD-(\d)"
        .- *.db: These two variables are use to relate the secuences that came 
        from  external databases. Using Dbxref.
    '''
    gff3 = []
    
    name = summary.query_name
    name = _transform_name_with_regexp(name, query_regex)
    query_length =  summary.query_length
    hits = summary.hits
    
    # Gbrowser needs a feature where the others features "hang on". It needs a
    # sequence with start and end to locate the other secuences. To match our 
    # Id schema, we have to trasform them with source database and regexps
    
    dbxref_query = _get_dbxref_name(query_db, name)   
    column_9  = "ID=%s,Name=%s%s " % (name, name, dbxref_query)
    feature_line = "%s\t.\t%s\t1\t%d\t.\t.\t.\t%s" % (name, feature_type,
                                                     query_length, column_9 )
    gff3.append(feature_line)
    
    
    for index, hit in enumerate(hits):
        hit_id = 'match_%s_%03d' % (name, index + 1)
        
        # here the order of the three functions is important !!
        hit_name = _transform_name_with_regexp(hit['name'], subject_regex)
        dbxref_subject = _get_dbxref_name(subject_db, hit_name)
        hit_name = _add_db_tag(hit_name, subject_db)
        
        column_9 = 'ID=%s;Name=%s%s' % (hit_id, hit_name, dbxref_subject)
        hit_line = '%s\t.\tmatch\t%d\t%d\t%s\t.\t.\t%s' % (name, hit['start'],
                                                      hit['end'],hit['expect'], 
                                                      column_9)
        gff3.append(hit_line)
        for hsp in hit['hsps']:
            column_9 = 'ID=%s;Name=%s;target=%s %d %d %s %s' % \
               (hit_id, hit_name, hit_name, hsp['subject_start'],
                 hsp['subject_end'], hsp['subject_strand'], dbxref_subject)
            hsp_line = '%s\t.\tmatch_part\t%d\t%d\t%s\t%s\t.\t%s' % \
              (name, hsp['query_start'], hsp['query_end'], hsp['expect'], 
               hsp['query_strand'], column_9)
            gff3.append(hsp_line)
    
    gff3.append('###')
    return "\n".join(gff3)
