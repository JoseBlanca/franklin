'''This module holds the code that allows to analyze the alignment search
result analysis.

It can deal with blasts, iprscan or ssaha2 results.
This results can be parsed, filtered and analyzed.
'''
from Bio.Blast import NCBIXML
from biolib.seqs import SeqWithQuality

from math import log10

class BlastParser(object):
    '''An iterator  blast parser that yields the blast results in a 
    multiblast file'''
    def __init__(self, fhand):
        'The init requires a file to be parser'
        fhand.seek(0, 0)
        self._blast_file  = fhand
        #we use the biopython parser
        self._blast_parse = NCBIXML.parse(fhand)

    def __iter__(self):
        'Part of the iterator protocol'
        return self

    @staticmethod
    def _create_result_structure(bio_result):
        'Given a BioPython blast result it returns our result structure'
        #the query name and definition
        #query = bio_result.query
        name  = bio_result.query_id
        definition = bio_result.query
        #try:
            #    name, definition = query.split(' ', 1)
            #except ValueError:
                #name = query
                #definition = None
        #length of query sequence
        length     = bio_result.query_letters 
        #now we can create the query sequence
        query = SeqWithQuality(name=name, description=definition,
                               length=length)
        
        #now we go for the hits (matches)
        matches = []
        for alignment in bio_result.alignments:
            #the subject sequence
            #subj_name = alignment.title.split()[0]
            name = alignment.accession
            definition = alignment.hit_def
            length = alignment.length
            subject = SeqWithQuality(name=name, description=definition,
                                     length=length)
            #the hsps (match parts)
            match_parts = []
            match_start, match_end = None, None 
            for hsp in alignment.hsps:
                expect         = hsp.expect
                subject_start  = hsp.sbjct_start
                subject_end    = hsp.sbjct_end
                query_start    = hsp.query_start
                query_end      = hsp.query_end
                hsp_length     = len(hsp.query)
                #We have to check the subject strand
                if subject_start < subject_end:
                    subject_strand = 1
                else:
                    subject_strand = -1
                    subject_start, subject_end = (subject_end,
                                                  subject_start)
                #Also the query strand
                if query_start < query_end:
                    query_strand = 1
                else:
                    query_strand = -1
                    query_start, query_end = query_end, query_start

                try:
                    similarity = hsp.positives*100.0/float(hsp_length)
                except TypeError:
                    similarity = None
                try:
                    identity = hsp.identities*100.0/float(hsp_length)
                except TypeError:
                    identity = None
                match_parts.append({
                    'subject_start'  : subject_start,
                    'subject_end'    : subject_end,
                    'subject_strand' : subject_strand,
                    'query_start'    : query_start,
                    'query_end'      : query_end,
                    'query_strand'   : query_strand,
                    'scores'         : {'similarity': similarity,
                                        'expect'    : expect,
                                        'identity'  : identity}
                    })
                # It takes the first loc and the last loc of the hsp to
                # determine hit start and end
                if match_start is None or query_start < match_start:
                    match_start = query_start
                if match_end is None or query_end > match_end:
                    match_end = query_end
                
            matches.append({
                'subject': subject,
                'start'  : match_start,   
                'end'    : match_end, 
                'scores' : {'expect':match_parts[0]['scores']['expect']},
                'match_parts' : match_parts})
        result = {'query'  : query,
                  'matches': matches}
        return result

    def next(self):
        'It returns the next blast result'
        bio_result = self._blast_parse.next()
        #now we have to change this biopython blast_result in our
        #structure
        our_result = self._create_result_structure(bio_result)
        return our_result

def _merge_overlaping_match_parts(match_parts):
    '''Given a list of match_parts it merges the ones that overlaps

       hsp 1  -------        ----->    -----------
       hsp 2       ------
       The similarity information will be lost
       The modified hsps will be in hsp_mod not in hsp
       It returns the list of new match_parts
    '''
    hsps = match_parts
    def cmp_subj_location(hsp_limit1, hsp_limit2):
        'It compares the subject locations'
        return hsp_limit1['subj'] - hsp_limit2['subj']
    
    hsp_limits = [] #all hsp starts and ends
    #we collect all start and ends
    for hsp in hsps:
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
    return merged_hsps

def _incompatible_length(match_parts, query, subject):
    '''It returns the incompatible length. The length of the region
    that should be alignned but it is not.

    It requires a list of match_parts, the query and the subject.
    '''
    first_hsp = match_parts[0]
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
    query_len   = len(query)
    subject_len = len(subject)

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

    last_hsp = match_parts[-1]
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
    if len(match_parts) == 1:
        return incompatible_len  #if there is only one hsp we're done
    to_add = 0

    #pylint: disable-msg=W0612
    for i, hsp in enumerate(match_parts[1:]):
        to_add = to_add + match_parts[i+1]['subject_start'] - \
                                        match_parts[i]['subject_end']
    incompatible_len += to_add
    return incompatible_len

class FilteredAlignmentResults(object):
    '''An iterator that yield the search results with its matches filtered
    
    It accepts a list of filters that will be applied to the alignment
    search results.
    '''
    def __init__(self, results, filters):
        '''It requires an alignment search result and a list of filters.
        
        The alignment search result is a dict with the query, matches,
        etc.
        The filters list can have severel filters defined by dicts.
        Each filter definition dict should have the key 'kind' as well
        as other keys that depend on the kind of filter.
        Allowed filters are:
            - best scores filters.
                It filters keeps the best match and the other matches
                equally good. Which is considered equally good is
                defined by the score_tolerance. All matches that
                obey:
                 (log10 best_match score - log10 match score) < tolerance
                are kept.
                    {'kind'           : 'best_scores',
                     'score_key'      : 'expect',
                     'min_score_value': 1e-4,
                     'score_tolerance': 10}
        '''
        self._results = results
        self._filters = filters

    def __iter__(self):
        'Part of the iteration protocol'
        return self

    def next(self):
        'It returns the next result filtered.'
        result = self._results.next()
        for filter_ in self._filters:
            self._filter_matches(result, filter_)
        return result

    @staticmethod
    def create_filter_best_score(parameters):
        'It returns a function that will filter matches'

        log_best_score = parameters['log_best_score']
        log_tolerance  = parameters['log_tolerance']
        score_key      = parameters['score_key']
        if 'min_score_value' in parameters:
            min_score  = parameters['min_score_value']
            max_score  = None
        else:
            min_score  = None
            max_score  = parameters['max_score_value']
        def filter_(match):
            '''It returns True or False depending on the match meeting
            the criteria'''
            #the score can be in the match itself or in the first
            #match_part
            if score_key in match['scores']:
                score = match['scores'][score_key]
            else:
                #the score is taken from the best hsp (the first one)
                score = match['match_parts'][0]['scores'][score_key]
            if max_score is not None and score == 0.0:
                result = True
            elif min_score is not None and score <= min_score:
                result = False
            elif max_score is not None and score >= max_score:
                result = False
            elif abs(log10(score) - log_best_score) < log_tolerance:
                result = True
            else:
                result = False
            return result
        return filter_

    @staticmethod
    def create_filter_min_score(parameters):
        'It returns a function that will filter matches'
        score_key      = parameters['score_key']
        if 'min_score_value' in parameters:
            min_score  = parameters['min_score_value']
            max_score  = None
        else:
            min_score  = None
            max_score  = parameters['max_score_value']
        def filter_(match):
            '''It returns True or False depending on the match meeting
            the criteria'''
            if score_key in match['scores']:
                score = match['scores'][score_key]
            else:
                #the score is taken from the best hsp (the first one)
                score = match['match_parts'][0]['scores'][score_key]
            if min_score is not None and score >= min_score:
                result = True
            elif max_score is not None and score <= max_score:
                result = True
            else:
                result = False
            return result
        return filter_

    @staticmethod
    def create_filter_min_length(parameters):
        'It returns a function that will filter matches'
        #the min length can be given in base pairs or as a percentage
        #of the query or the subject
        kind, query, min_length = None, None, None
        if 'min_length_bp' in parameters:
            min_length = parameters['min_length_bp']
            kind = 'bp'
        elif 'min_length_query_%' in parameters:
            min_length = parameters['min_length_query_%']
            kind  = 'query'
            query = parameters['query']
        elif 'min_length_subject_%' in parameters:
            min_length = parameters['min_length_subject_%']
            kind = 'subject'
        else:
            raise ValueError('Filter poorly defined, missing parameters')

        def filter_(match):
            '''It returns True or False depending on the match meeting
            the criteria'''
            #how to calculate the match length depends on the
            #kind of filtering we're doing: base pairs, percentage
            #on the query or on the subject
            match_length = match['end'] - match['start'] + 1
            if kind == 'bp':
                match_length = match_length
            elif kind == 'query':
                match_length = (match_length / float(len(query))) * 100.0
            elif kind == 'subject':
                subject = match['subject']
                match_length = \
                             (match_length / float(len(subject))) * 100.0
            if match_length >= min_length:
                result = True
            else:
                result = False
            return result
        return filter_

    @staticmethod
    def create_filter_compatibility(parameters):
        '''It returns a function that will filter matches
        It select only the ones that are compatible in more than
        min_compatibility. Compatible is the fragment that is aligned in
        a match_part. Also it shouldn't be more than max_incompatibility.
        Incompatibility means that they should be alignmed in that region,
        but they aren't. For instance:
            -----------------------------------
                                |||||||||||||||
                -------------------------------
                incompatible    compatible
                  region          region
        Both compatibility and incompatibility will be applied to both
        query and subject and are number of residues.
        '''
        #the min length can be given in base pairs or as a percentage
        #of the query or the subject
        min_compat   = parameters['min_compatibility']
        max_incompat = parameters['max_incompatibility']
        min_simil    = parameters['min_similarity']
        query        = parameters['query']

        def filter_(match):
            '''It returns True or False depending on the match meeting
            the criteria'''
            #filter match_parts under similarity
            match_parts = []
            for match_part in match['match_parts']:
                if match_part['scores']['similarity'] >= min_simil:
                    match_parts.append(match_part)
            if not match_parts:
                return False
            #merge_overlaping_match_parts
            match_parts = _merge_overlaping_match_parts(match_parts)

            #are the remaining match_parts lengthier than the min?
            length = 0
            for match_part in match_parts:
                length = length + match_part['subject_end'] - \
                                          match_part['subject_start'] + 1
            if length < min_compat:
                return False
            #filter matchs above incompatibility threshold
            subject = match['subject']
            length = _incompatible_length(match_parts, query, subject)
            if length > max_incompat:
                return False
            return True
        return filter_

    def _filter_matches(self, result, filter_):
        'Given a filter dict and a result it filters the matches'
        filter_functions_factories = {
            'best_scores'  : self.create_filter_best_score,
            'min_scores'   : self.create_filter_min_score,
            'min_length'   : self.create_filter_min_length,
            'compatibility': self.create_filter_compatibility,
        }

        kind = filter_['kind']

        #some filters need extra data
        if kind == 'best_scores':
            #the log10 for the best score
            #the best score should be in the first hit
            score_key = filter_['score_key']
            best_score = result['matches'][0]['scores'][score_key]
            if best_score == 0.0:
                log_best_score = 0.0
            else:
                log_best_score = log10(best_score)
            filter_['log_best_score'] = log_best_score
            filter_['log_tolerance']  = log10(filter_['score_tolerance'])
        elif kind == 'min_length' and 'min_length_query_%' in filter_:
            filter_['query'] = result['query']
        elif kind == 'compatibility':
            filter_['query'] = result['query']

        filter_ = filter_functions_factories[kind](filter_)
        #pylint: disable-msg=W0141
        result['matches'] = list(filter(filter_, result['matches']))

