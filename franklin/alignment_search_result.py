'''This module holds the code that allows to analyze the alignment search
result analysis.

It can deal with blasts, iprscan or ssaha2 results.
This results can be parsed, filtered and analyzed.

This module revolves around a memory structure that represents a blast or
an iprscan result. The schema of this structure is:
result = {'query':the_query_sequence,
          'matches': [a_list_of_matches(hits in the blast terminology)]
         }
The sequence can have: name, description, annotations={'database':some db} and
len(sequence).
Every match is a dict.
match  = {'subject':the subject sequence
          'start'  :match start position in bp in query
          'end'    :match end position in bp in query
          'subject_start' : match start position in bp in subject
          'subject_end'    :match end position in bp in subject
          'scores' :a dict with the scores
          'match_parts': [a list of match_parts(hsps in the blast lingo)]
          'evidences'  : [a list of tuples for the iprscan]
         }
All the scores are holded in a dict
scores  = {'key1': value1, 'key2':value2}
For instance the keys could be expect, similarity and identity for the blast

match_part is a dict:
    match_part = {'query_start'    : the query start in the alignment in bp
                  'query_end'      : the query end in the alignment in bp
                  'query_strand'   : 1 or -1
                  'subject_start'  : the subject start in the alignment in bp
                  'subject_end'    : the subject end in the alignment in bp
                  'subject_strand' : 1 or -1
                  'scores'         :a dict with the scores
            }
Iprscan has several evidences generated by different programs and databases
for every match. Every evidence is similar to a match.
'''

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

from Bio.Blast import NCBIXML
from franklin.seq.seqs import SeqWithQuality
from Bio.Seq import UnknownSeq
from math import log10

class BlastParser(object):
    '''An iterator  blast parser that yields the blast results in a
    multiblast file'''
    def __init__(self, fhand):
        'The init requires a file to be parser'
        fhand.seek(0, 0)
        self._blast_file = fhand
        blast_version, plus = self._get_blast_version()
        self._blast_file.seek(0, 0)
        
        if ((blast_version and plus) or 
                                (blast_version and blast_version >= '2.2.21')):
            self.use_query_def_as_accession = True
            self.use_subject_def_as_accession = True
        
        else:
            self.use_query_def_as_accession = True
            self.use_subject_def_as_accession = False

        #we use the biopython parser
        #if there are no results we put None in our blast_parse results
        self._blast_parse = None
        if fhand.read(1) == '<':
            fhand.seek(0)
            self._blast_parse = NCBIXML.parse(fhand)

    def __iter__(self):
        'Part of the iterator protocol'
        return self


    def _create_result_structure(self, bio_result):
        'Given a BioPython blast result it returns our result structure'
        #the query name and definition
        definition = bio_result.query
        if self.use_query_def_as_accession:
            items = definition.split(' ', 1)
            name = items[0]
            if len(items) > 1:
                definition = items[1]
            else:
                definition = None
        else:
            name = bio_result.query_id
            definition = definition
        if definition is None:
            definition = "<unknown description>"
        #length of query sequence
        length = bio_result.query_letters
        #now we can create the query sequence
        query = SeqWithQuality(name=name, description=definition,
                               seq=UnknownSeq(length=length))

        #now we go for the hits (matches)
        matches = []
        for alignment in bio_result.alignments:
            #the subject sequence
            if self.use_subject_def_as_accession:
                items = alignment.hit_def.split(' ', 1)
                name = items[0]
                if len(items) > 1:
                    definition = items[1]
                else:
                    definition = None
            else:
                name = alignment.accession
                definition = alignment.hit_def

            if definition is None:
                definition = "<unknown description>"

            length = alignment.length
            subject = SeqWithQuality(name=name, description=definition,
                                        seq=UnknownSeq(length=length))

            #the hsps (match parts)
            match_parts = []
            match_start, match_end = None, None
            match_subject_start, match_subject_end = None, None
            for hsp in alignment.hsps:
                expect = hsp.expect
                subject_start = hsp.sbjct_start
                subject_end = hsp.sbjct_end
                query_start = hsp.query_start
                query_end = hsp.query_end
                hsp_length = len(hsp.query)
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
                    similarity = hsp.positives * 100.0 / float(hsp_length)
                except TypeError:
                    similarity = None
                try:
                    identity = hsp.identities * 100.0 / float(hsp_length)
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
                if match_subject_start is None or subject_start < match_subject_start:
                    match_subject_start = subject_start
                if match_subject_end is None or subject_end > match_subject_end:
                    match_subject_end = subject_end
            matches.append({
                'subject': subject,
                'start'  : match_start,
                'end'    : match_end,
                'subject_start': match_subject_start,
                'subject_end'  : match_subject_end,
                'scores' : {'expect':match_parts[0]['scores']['expect']},
                'match_parts' : match_parts})
        result = {'query'  : query,
                  'matches': matches}
        return result
    def _get_blast_version(self):
        'It gets blast parser version'
        for line in self._blast_file:
            line = line.strip()
            if line.startswith('<BlastOutput_version>'):
                version = line.split('>')[1].split('<')[0].split()[1]
                break
        plus = False
        if '+' in version:
            plus = True
            version = version[:-1] 
        return version, plus
        

    def next(self):
        'It returns the next blast result'
        if self._blast_parse is None:
            raise StopIteration
        else:
            bio_result = self._blast_parse.next()
            #now we have to change this biopython blast_result in our
            #structure
            our_result = self._create_result_structure(bio_result)
            return our_result

class ExonerateParser(object):
    '''Exonerate parser, it is a iterator that yields the result for each
    query separated'''

    def __init__(self, fhand):
        'The init requires a file to be parser'
        self._fhand = fhand
        self._exonerate_results = self._results_query_from_exonerate()
    def __iter__(self):
        'Part of the iterator protocol'
        return self

    def _results_query_from_exonerate(self):
        '''It takes the exonerate cigar output file and yields the result for
        each query. The result is a list of match_parts '''
        self._fhand.seek(0, 0)
        cigar_dict = {}
        for line in  self._fhand:
            if not line.startswith('cigar_like:'):
                continue
            items = line.split(':', 1)[1].strip().split()
            query_id = items[0]
            if query_id not in cigar_dict:
                cigar_dict[query_id] = []
            cigar_dict[query_id].append(items)
        for query_id, values in cigar_dict.items():
            yield values

    @staticmethod
    def _create_structure_result(query_result):
        '''It creates the result dictionary structure giving a list of
        match_parts of a query_id '''
        #TODO add to the match the match subject start and end
        struct_dict = {}
        query_name = query_result[0][0]
        query_length = int(query_result[0][9])

        query = SeqWithQuality(name=query_name,
                                      seq=UnknownSeq(length=query_length))
        struct_dict['query'] = query
        struct_dict['matches'] = []
        for match_part_ in query_result:
            (query_name, query_start, query_end, query_strand, subject_name,
            subject_start, subject_end, subject_strand, score, query_length,
            subject_length, similarity) = match_part_
            query_start = int(query_start)
            #they number the positions between symbols
            # A C G T
            #0 1 2 3 4
            #Hence the subsequence "CG" would have start=1, end=3, and length=2
            #but we would say start=1 and end=2
            query_end = int(query_end) - 1
            subject_start = int(subject_start)
            subject_end = int(subject_end) - 1
            query_strand = _strand_transform(query_strand)
            subject_strand = _strand_transform(subject_strand)
            score = int(score)
            similarity = float(similarity)
            # For each line , It creates a match part dict
            match_part = {}
            match_part['query_start'] = query_start
            match_part['query_end'] = query_end
            match_part['query_strand'] = query_strand
            match_part['subject_start'] = subject_start
            match_part['subject_end'] = subject_end
            match_part['subject_strand'] = subject_strand
            match_part['scores'] = {'score':score,
                                           'similarity':similarity}

            # Check if the match is already added to the struct. A match is
            # defined by a list of part matches between a query and a subject
            match_num = _match_num_if_exists_in_struc(subject_name, struct_dict)
            if match_num is not None:
                match = struct_dict['matches'][match_num]
                if match['start'] > query_start:
                    match['start'] = query_start
                if match['end'] < query_end:
                    match['end'] = query_end
                if match['scores']['score'] < score:
                    match['scores']['score'] = score
                match['match_parts'].append(match_part)
            else:
                match = {}
                match['subject'] = SeqWithQuality(name=subject_name,
                                     seq=UnknownSeq(length=int(subject_length)))
                match['start'] = query_start
                match['end'] = query_end
                match['scores'] = {'score':score}
                match['match_parts'] = []
                match['match_parts'].append(match_part)
                struct_dict['matches'].append(match)
        return struct_dict

    def next(self):
        '''It return the next exonerate hit'''
        query_result = self._exonerate_results.next()
        return self._create_structure_result(query_result)

def _strand_transform(strand):
    '''It transfrom the +/- strand simbols in our user case 1/-1 caracteres '''
    if strand == '-':
        return - 1
    elif strand == '+':
        return 1

def _match_num_if_exists_in_struc(subject_name, struct_dict):
    'It returns the match number of the list of matches that is about subject'
    for i, match in enumerate(struct_dict['matches']):
        if subject_name == match['subject'].name:
            return i
    return None

def incompatibility_score(match, query, subject):
    '''It calculates and returns the incompatibility score

    The incompatible region is the region that should be aligned, but it is
    not.

       query     -----------------
                   |||||
       subject   ----------------
                 ++     +++++ <-incompabible regions
    The score is given as a percentage relative to the shortest sequence
    beteween the query and the subject.
    '''
    incomp = _compatible_incompatible_length(match, query)[1]
    #we calculate a percentage dividing by the shortest seq
    #between query an subj
    min_len = min(len(query), len(subject))
    incomp_bak = incomp
    incomp = float(incomp) / float(min_len) * 100.0
    if incomp < 0 or incomp > 100:
        print 'incomp', incomp, 'match', match
        print 'min_len, incomp, %', min_len, incomp_bak, incomp
        msg = 'Bad calculation for incompatible region'
        raise RuntimeError(msg)
    return incomp

def get_match_score(match, score_key, query=None, subject=None):
    '''Given a match it returns its score.

    It tries to get the score from the match, if it's not there it goes for
    the first match_part.
    It can also be a derived score like the incompatibility. All derived scores
    begin with d_
    '''
    derived_scores = {'d_incompatibility': incompatibility_score}
    #the score can be in the match itself or in the first
    #match_part
    if score_key in match['scores']:
        score = match['scores'][score_key]
    elif score_key in derived_scores:
        score = derived_scores[score_key](match, query, subject)
    else:
        #the score is taken from the best hsp (the first one)
        score = match['match_parts'][0]['scores'][score_key]
    return score

def _merge_overlaping_match_parts(match_parts, min_similarity=None):
    '''Given a list of match_parts it merges the ones that overlaps

       hsp 1  -------        ----->    -----------
       hsp 2       ------
       The similarity information will be lost
       The modified hsps will be in hsp_mod not in hsp
       It returns the list of new match_parts
    '''
    #DO NOT USE THIS FUNCTION, LOTS OF PATOLOGICAL CASES!!!!!!
    #we don't use the merging of all match_parts because is very difficult
    #to calculate the correct one
    #what happens in a case like:
    #    subj 1-----------78           478-------------534
    #   query 1-----------78           31--------------78
    hsps = match_parts
    #hsp
    hsp0 = hsps[0]
    subject_strand = hsp0['subject_strand']
    query_strand = hsp0['query_strand']
    hsp0_qe = hsp0['query_end']
    hsp0_se = hsp0['subject_end']

    #new hsps
    filtered_hsps = [hsp0]
    for hsp in hsps[1:]:
        if min_similarity and hsp['scores']['similarity'] < min_similarity:
            continue
        if (hsp['subject_strand'] != subject_strand or
            hsp['query_strand'] != query_strand):
            #we'll have into account only the matches with the same
            #orientation as the first hsp
            continue
        #what happens in a case like:
        #    subj 1-----------78           478-------------534
        #   query 1-----------78           31--------------78
        #so we only use the hsps that are in the same diagonal as the first
        if query_strand == 1:
            hsp_qs = hsp['query_start']
            hsp_ss = hsp['subject_start']
            if hsp_ss - hsp0_se != 0:
                slope = float(hsp_qs - hsp0_qe) / float(hsp_ss - hsp0_se)
            else:
                slope = float(hsp['query_end'] - hsp0_qe) / \
                        float(hsp['subject_end'] - hsp0_se)
        else:
            hsp_ss = hsp['subject_start']
            if hsp0_se - hsp_ss != 0:
                slope = float(hsp0['query_start'] - hsp['query_end']) / \
                        float(hsp0_se - hsp_ss)
            else:
                slope = float(hsp0['query_start'] - hsp['query_start']) / \
                        float(hsp['subject_end'] - hsp0_se)
        #slope limits
        slope_min = 0.95
        slope_max = 1.05
        if query_strand != subject_strand:
            slope_min, slope_max = -slope_max, -slope_min
        if slope >= slope_min and slope <= slope_max:
            filtered_hsps.append(hsp)
    hsps = filtered_hsps
    #print 'hsps', len(hsps)

    #we collect all start and ends
    hsp_limits = [] #all hsp starts and ends
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

    #now we sort the hsp limits according their query location
    def cmp_query_location(hsp_limit1, hsp_limit2):
        'It compares the query locations'
        return hsp_limit1['query'] - hsp_limit2['query']
    hsp_limits.sort(cmp_query_location)

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
                query_strand = None
                merged_hsps.append({
                    'expect' : None,
                    'subject_start' : subj_start,
                    'subject_end' : subj_end,
                    'subject_strand':hsp0['subject_strand'],
                    'query_start' : query_start,
                    'query_end' : query_end,
                    'query_strand':hsp0['query_strand'],
                    'similarity' : None
                    }
                )
    return merged_hsps

def _compatible_incompatible_length(match, query, min_similarity=None):
    '''It returns the compabible and incompatible length in a match

        xxxxxxx      xxxxxxxxxx      xxx incompatible
    ------------------------------------
               ||||||          ||||||
        ------------------------------------------
               ******          ******    compatible
    '''
    match_parts = _merge_overlaping_match_parts(match['match_parts'],
                                                min_similarity)
    #     first_incomp
    #     xxxxxxx      xxxxxxxxxx      xxx incompatible
    #------------------------------------
    #           ||||||          ||||||
    #    ------------------------------------------
    #           ******          ******    compatible
    first_matchp = match_parts[0]
    if first_matchp['query_strand'] == 1:
        first_incomp = min(first_matchp['query_start'],
                           first_matchp['subject_start'])
    else:
        first_incomp = min(first_matchp['query_start'],
                           len(match['subject']) - \
                                             first_matchp['subject_end'] + 1)

    #print '1_match', first_matchp
    last_matchp = match_parts[-1]
    query_overlap = len(query) - last_matchp['query_end'] - 1
    if first_matchp['query_strand'] == 1:
        subject_overlap = len(match['subject']) - last_matchp['subject_end'] \
                                                                           - 1
    else:
        subject_overlap = last_matchp['subject_start']
    last_incomp = min(query_overlap, subject_overlap)

    #the compatible length for all match_parts
    comp_length = 0
    for matchp in match_parts:
        comp_length += matchp['query_end'] - matchp['query_start'] + 1

    #which is the incompatible regions between consecutive match parts?
    match_incomp = last_matchp['query_end'] - first_matchp['query_start'] + \
                                                               1 - comp_length
    #print '1, last, match', first_incomp, last_incomp, match_incomp
    incomp = first_incomp + last_incomp + match_incomp
    return comp_length, incomp

def get_alignment_parser(kind):
    '''It returns a parser depending of the aligner kind '''
    parsers = {'blast':BlastParser, 'exonerate':ExonerateParser,
               'blast+':BlastParser}
    parser = parsers[kind]
    return parser

class FilteredAlignmentResults(object):
    '''An iterator that yield the search results with its matches filtered

    It accepts a list of filters that will be applied to the alignment
    search results.
    '''
    def __init__(self, results, match_filters=None, result_filters=None):
        '''It requires an alignment search result and a list of filters.

        The alignment search result is a dict with the query, matches,
        etc.
        The match_filters list can have several filters defined by dicts.
        Each filter definition dict should have the key 'kind' as well
        as other keys that depend on the kind of filter.
        Allowed match_filters are:
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

        The results filters list is a list of filters that can filter the
        results not the matches.
            Example:It can filter results by number of matches
        '''
        self._results = results
        if match_filters is None:
            match_filters = []
        self._match_filters = match_filters
        if result_filters is None:
            result_filters = []
        #for the results we create a default filter, we don't want results with
        #no matches
        result_filters.append({'kind': 'min_num_matches', 'value': 1})
        self._result_filters = result_filters

    def __iter__(self):
        'Part of the iteration protocol'
        return self

    def next(self):
        'It returns the next result filtered.'
        while True:
            result = self._results.next()
            #apply the match filters to the results
            for filter_ in self._match_filters:
                self._filter_matches(result, filter_)

            #is the result worth? we apply the result filters to find out.
            result_filters = {
                'max_num_matches':     self._max_num_matches_result_filter,
                'min_num_matches':     self._min_num_matches_result_filter,
                'max_num_match_parts': self._max_num_match_parts_result_filter}

            for result_filter in self._result_filters:
                filter_funct = result_filters[result_filter['kind']]
                if not filter_funct(result, result_filter):
                    result = None
                    break

            #is there a result to return or it has been filtered out?
            if result is not None:
                return result

    @staticmethod
    def _max_num_matches_result_filter(result, parameters):
        'It returns true if the given result less or equal matches'
        if len(result['matches']) <= parameters['value']:
            return True
        else:
            return False

    @staticmethod
    def _max_num_match_parts_result_filter(result, parameters):
        'It returns true if the given result less or equal matches'
        num_match_parts = 0
        for match in result['matches']:
            num_match_parts += len(match['match_parts'])
        if num_match_parts <= parameters['value']:
            return True
        else:
            return False

    @staticmethod
    def _min_num_matches_result_filter(result, parameters):
        'It returns true if the given result less or equal matches'
        if len(result['matches']) >= parameters['value']:
            return True
        else:
            return False

    @staticmethod
    def _create_filter_best_score(parameters):
        'It returns a function that will filter matches'

        log_best_score = parameters['log_best_score']
        log_tolerance = parameters['log_tolerance']
        score_key = parameters['score_key']
        if 'min_score_value' in parameters:
            min_score = parameters['min_score_value']
            max_score = None
        else:
            min_score = None
            max_score = parameters['max_score_value']
        def filter_(match):
            '''It returns True or False depending on the match meeting
            the criteria'''
            score = get_match_score(match, score_key)
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
    def _create_filter_min_score(parameters):
        'It returns a function that will filter matches'
        score_key = parameters['score_key']
        if 'min_score_value' in parameters:
            min_score = parameters['min_score_value']
            max_score = None
        else:
            min_score = None
            max_score = parameters['max_score_value']
        def filter_(match):
            '''It returns True or False depending on the match meeting
            the criteria'''
            score = get_match_score(match, score_key)
            if min_score is not None and score >= min_score:
                result = True
            elif max_score is not None and score <= max_score:
                result = True
            else:
                result = False
            return result
        return filter_

    @staticmethod
    def _create_filter_min_length(parameters):
        'It returns a function that will filter matches'
        #the min length can be given in base pairs or as a percentage
        #of the query or the subject
        kind, query, min_length = None, None, None
        if 'min_length_bp' in parameters:
            min_length = parameters['min_length_bp']
            kind = 'bp'
        elif 'min_length_query_%' in parameters:
            min_length = parameters['min_length_query_%']
            kind = 'query'
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
    def _create_filter_compatibility(parameters):
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
        min_compat = parameters['min_compatibility']
        max_incompat = parameters['max_incompatibility']
        min_simil = parameters['min_similarity']
        query = parameters['query']

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
            #which spans are aligned and which not?
            compat, incompat = _compatible_incompatible_length(match, query,
                                                               min_simil)
            if compat < min_compat:
                return False
            if incompat > max_incompat:
                return False
            return True
        return filter_

    def _filter_matches(self, result, filter_):
        'Given a filter dict and a result it filters the matches'
        filter_functions_factories = {
            'best_scores'  : self._create_filter_best_score,
            'min_scores'   : self._create_filter_min_score,
            'min_length'   : self._create_filter_min_length,
            'compatibility': self._create_filter_compatibility,
        }

        kind = filter_['kind']

        #some filters need extra data
        if kind == 'best_scores':
            score_key = filter_['score_key']
            #if the are no matches, there's nothig to filter
            if not result['matches']:
                return
            #the log10 for the best score
            #the best score should be in the first hit
            best_score = result['matches'][0]['scores'][score_key]
            if best_score == 0.0:
                log_best_score = 0.0
            else:
                log_best_score = log10(best_score)
            filter_['log_best_score'] = log_best_score
            filter_['log_tolerance'] = log10(filter_['score_tolerance'])
        elif kind == 'min_length' and 'min_length_query_%' in filter_:
            filter_['query'] = result['query']
        elif kind == 'compatibility':
            filter_['query'] = result['query']

        filter_ = filter_functions_factories[kind](filter_)
        #pylint: disable-msg=W0141
        result['matches'] = list(filter(filter_, result['matches']))

def get_match_scores(match, score_keys, query, subject):
    '''It returns the scores for one match.

    scores should be a list and it will return a list of scores.
    '''
    scores_res = []
    for score_key in score_keys:
        score = get_match_score(match, score_key, query, subject)
        scores_res.append(score)
    return scores_res

def alignment_results_scores(results, scores, filter_same_query_subject=True):
    '''It returns the list of scores for all results.

    For instance, for a blast a generator with all e-values can be generated.
    By default, the results with the same query and subject will be filtered
    out.
    The scores can be a single one or a list of them.
    '''
    #for each score we want a list to gather the results
    score_res = []
    for score in scores:
        score_res.append([])
    for result in results:
        query = result['query']
        for match in result['matches']:
            subject = match['subject']
            if (filter_same_query_subject and query is not None and subject is
                not None and query.name == subject.name):
                continue
            #all the scores for this match
            score_values = get_match_scores(match, scores, query, subject)
            #we append each score to the corresponding result list
            for index, value in enumerate(score_values):
                score_res[index].append(value)
    if len(score_res) == 1:
        return score_res[0]
    else:
        return score_res

'''
A graphical overview of a blast result can be done counting the hits
with a certain level of similarity. We can represent a distribution
with similarity in the x axe and the number of hits for each
similarity in the y axe. It would be something like.

 n |
 u |
 m |
   |       x
 h |      x x
 i |     x   x     x
 t |    x      x  x x
 s | xxx        xx   x
    ----------------------
         % similarity

Looking at this distribution we get an idea about the amount of
similarity found between the sequences used in the blast.

Another posible measure between a query and a subject is the
region that should be aligned, but it is not, the incompatible
region.

   query         -----------------
                   |||||
   subject   ----------------
                 ++     +++++ <-incompabible regions

For every query-subject pair we can calculate the similarity and the
incompatible region. We can also draw a distribution with both
measures. (a 3-d graph viewed from the top with different colors
for the different amounts of hits).

 % |
 s |
 i |       ..
 m |
 i |             ..
 l |           ..
 a |     x
 r |      xx         ..
 i |
 t | xx       xx
 y | xx       xx
   -----------------------
    % incompatibility

'''
