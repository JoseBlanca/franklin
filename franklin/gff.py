'''This module deals with the IO of GFF files.

Created on 26/10/2009
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

from franklin.seq.readers import seqs_in_file
from franklin.seq.writers import write_seqs_in_file
ESCAPES = {';': '%3B',
           '=': '%3D',
           '&': '%26',
           ' ': '%20',
           '%': '%25',
           ',': '%2C'}

_DEFAULT_READ_VERSION = '2'
_DEFAULT_WRITE_VERSION = '3'
METADATA = 'metadata'
COMMENT = 'comment'
FEATURE = 'feature'
FASTA = 'fasta'
_ATTRIBUTE_KEYVAL_SEPARATORS = {'2': ' ', '3':'='}

def _deescape(string):
    'It replaces the escaped sequences with the real characters'

    for character, escape in ESCAPES.items():
        string = string.replace(escape, character)
    return string

class GffFile(object):
    'A GFF object for reading and writing'
    def __init__(self, fpath, mode='r'):
        '''It inits the class,

        Accepted modes are r and w.
        '''
        self._mode = mode
        self._fpath = fpath

        if mode not in ('r', 'w'):
            msg = 'Mode should be r or w. Mode not supported: %s' % mode
            raise ValueError(msg)
        if mode == 'w':
            self._fhand = open(fpath, mode)
        else:
            self._fhand = None

        #if we're reading, is the version on the first line?
        self._version = None
        if mode == 'r':
            self._read_version()
        self._feature_ids = {}

    def _read_version(self):
        'It reads the GFF version'
        line = open(self._fpath, self._mode).readline()

        if line.startswith('##gff-version'):
            self._version = line.strip().split()[1]
        else:
            self._version = _DEFAULT_READ_VERSION

    def _get_version(self):
        'It returns the GFF version'
        return self._version
    def _set_version(self, version):
        'It allows to change the version if no line has been written'
        # we should be able to set version when we are reading.
        # this is for gff version 3 without specific header ##gff-version 3
        # like in version 2011-03-11
        if self._fhand is None:
            pass
        elif self._fhand.tell() != 0:
            msg = 'It is not possible to set the version after writing a line'
            raise RuntimeError(msg)
        self._version = version
    version = property(_get_version, _set_version)

    def _create_feature_id(self, feature_name):
        'It returns a new unique feature_id based on the name'
        feature_ids = self._feature_ids
        if feature_name not in feature_ids:
            feature_id = feature_name
            feature_ids[feature_name] = 1
        else:
            feature_id = '%s_%i' % (feature_name, feature_ids[feature_name])

        feature_ids[feature_name] += 1
        return feature_id

    def _create_feature(self, line):
        'It creates a feature from a line of a gff'
        feature = {}
        try:
            seqid, source, type_, start, end, score, strand, phase, annots = \
                                                            line.split("\t" , 8)
        except ValueError:
            msg = 'Malformed feature: %s' % line
            raise ValueError(msg)
        start = int(start)
        end   = int(end)

        attributes = {}
        separator = _ATTRIBUTE_KEYVAL_SEPARATORS[self.version]

        feature_id = None
        for attribute in annots.rstrip(';').split(';'):
            attribute = attribute.strip(' ')
            try:
                key, value = attribute.split(separator, 1)
            except ValueError:
                msg = 'Malformed attribute: %s' % attribute
                raise ValueError(msg)

            value = value.strip('"')
            key = _deescape(key)
            value = _deescape(value)
            attributes[key] = value
            if key == 'Name':
                feature['name'] = value
            if key == 'ID':
                feature_id = value

        if not feature_id:
            feature_id = self._create_feature_id(feature['name'])

        feature['id'] = _deescape(feature_id)
        feature['seqid']   = _deescape(seqid)
        feature['source']  = _deescape(source)
        feature['type']    = type_
        feature['start']   = start
        feature['end']     = end
        feature['score']   = score
        feature['strand']  = strand
        feature['phase']   = phase
        feature['attributes'] = attributes

        return feature

    def _get_items(self):
        'It yields the items in the GFF file'

        fhand = open(self._fpath, self._mode)

        for line in fhand:
            line = line.strip()
            if not line:
                continue
            if line.startswith('##gff-version'):
                continue #this has been taken into account
            elif line.startswith('##FASTA'):
                #in the next line we're assuming that the fasta reading is
                #sequential and that there is no seek
                items = seqs_in_file(fhand, format='fasta')
                yield FASTA, items
                break
            elif line.startswith('##'):
                item = line[2:]
                kind = METADATA
            elif line.startswith('#'):
                item = line[1:]
                kind = COMMENT
            else:
                item = self._create_feature(line)
                kind = FEATURE
            if item is not None:
                yield kind, item

    items = property(_get_items)

    def _get_features(self):
        'It yields the features in the GFF file'
        for kind, item in self.items:
            if kind == FEATURE:
                yield item
    features = property(_get_features)

    def write(self, item):
        '''It writes a line.

        The item should be a tuple with the kind and the information about the
        feature
        '''
        if item is None:
            return
        kind, item = item
        if self._fhand.tell() == 0:
            if not self.version:
                self.version = _DEFAULT_WRITE_VERSION
            self._fhand.write('##gff-version %s\n' % self.version)
        if kind == METADATA:
            self._fhand.write('##' + item + '\n')
        elif kind == COMMENT:
            self._fhand.write('#' + item + '\n')
        elif kind == FEATURE:
            feature_line = self._feature_to_str(item) + '\n'
            self._fhand.write(feature_line)
        elif kind == FASTA:
            self._fhand.write('##FASTA\n')
            write_seqs_in_file(item, self._fhand, format='fasta')
        self._fhand.flush()

    @staticmethod
    def _escape(string, escape_coma=True):
        'It returns an escaped string'
        #codes taken from:
        #http://www.blooberry.com/indexdot/html/
        #topics/urlencoding.htm?state=urlenc&origval=%5Ct&enc=on
        escapes = ESCAPES.copy()
        escape_strings = escapes.values()
        if not escape_coma:
            del escapes[',']

        new_string = []
        for index, char_ in enumerate(string):
            #it might be already escaped
            if char_ == '%':
                is_escaped = False
                try:
                    full_escape = string[index:index + 3]
                    if full_escape in escape_strings:
                        is_escaped = True
                except IndexError:
                    new_string.append(char_)
                if not is_escaped:
                    char_ = escapes[char_]
            elif char_ in escapes:
                char_ = escapes[char_]
            new_string.append(char_)
        return ''.join(new_string)

    def _feature_to_str(self, feature):
        'Given a feature dict it returns a gff feature line'
        feature_fields = []

        feature_fields.append(self._escape(feature['seqid']))

        if 'source' in feature:
            feature_fields.append(self._escape(feature['source']))
        else:
            feature_fields.append('.')

        feature_fields.append(feature['type'])
        feature_fields.append('%d' % int(feature['start']))
        feature_fields.append('%d' % int(feature['end']))

        for tag, default in (('score', '.'), ('strand', '.'), ('phase', '.')):
            if tag in feature:
                feature_fields.append(feature[tag])
            else:
                feature_fields.append(default)

        # attributes
        if 'attributes' in feature:
            attributes = feature['attributes']
        else:
            attributes = {}

        if 'name' in feature:
            feat_name = feature['name']
        elif 'name' in attributes:
            feat_name = attributes['name']
        elif 'Name' in attributes:
            feat_name = attributes['Name']
        attributes['Name'] = feat_name

        if 'id' in feature:
            feature_id = feature['id']
        elif 'ID' in attributes:
            feature_id = attributes['ID']
        else:
            feature_id = self._create_feature_id(attributes['Name'])

        attributes['ID'] = feature_id

        attribute_list = []
        for attr_key, attr_value in attributes.items():
            if isinstance(attr_value, list):
                attr_value = ','.join(attr_value)
            if attr_value is None:
                continue
            attr_key = self._escape(attr_key)
            attr_value = self._escape(attr_value, escape_coma=False)
            attribute_list.append('%s=%s' % (attr_key, attr_value))

        feature_fields.append(';'.join(attribute_list))

        return '\t'.join(feature_fields)

    def flush(self):
        'It flushes the file'
        self._fhand.flush()

def write_gff(out_fpath, items):
    '''It writes a GFF file.

    It requires a list of items.
    Every item is a tuple with a kind and a content.
    Allowed kinds are: metadata, comment, feature and fasta
    A feature is a dict with the following keys: seqid, source, kind, start,
    end, score, strand, phase, id, name, alias, target, gap, derives_from, note,
    dbxref, ontology_term, parents, and children.
    '''
    writer = GffFile(fpath=out_fpath, mode='w')

    for kind, item in items:
        writer.write((kind, item))
    writer.flush()

def create_feature_type_filter(types):
    'it creates a filter that filters by type of the feature'

    def feature_type_filter(item):
        'The real filter'
        if item[0] != FEATURE:
            return item
        kind = item[1]['type']
        if kind in types:
            return True
        else:
            return False
    return feature_type_filter

def create_mapper_add_id_as_name():
    'It adds the id as name if name not in gff'
    def id_to_name_mapper(item):
        'It adds the name using id if not name'
        if item[0] != FEATURE:
            return item
        feature = item[1]
        if 'name' not in feature:
            feature['name'] = feature['id']
        return (FEATURE, feature)
    return id_to_name_mapper

def create_description_adder(descriptions):
    '''it adds the description from a dict.

    It adds it on the Note field of attributes'''
    def add_annot_mapper(item):
        'It adds descriptions to the feature'
        if item[0] != FEATURE:
            return item
        feature = item[1]
        attr_key = 'Note'
        name = feature.get('name', None)
        if not name:
            return feature
        description = descriptions.get(name, None)
        if description:
            attr =_add_dbxrefs_to_dbxref(feature['attributes'].get(attr_key, ''),
                                         [description])
            feature['attributes'][attr_key] = attr
        return (FEATURE, feature)
    return add_annot_mapper

def create_go_annot_adder(go_terms):
    '''It creates a mapper that adds gene ontology terms.

    The GO terms should be provided by a dict with the feature names as keys
    and a list of GO terms as values
    '''
    def go_annot_mapper(item):
        'It adds the GO terms to the feature'
        if item[0] != FEATURE:
            return item
        feature = item[1]
        attr_key = 'Ontology_term'
        name = feature.get('name', None)
        if not name:
            return (FEATURE, feature)
        gos = go_terms.get(name, None)
        if gos:
            attr =_add_dbxrefs_to_dbxref(feature['attributes'].get(attr_key, ''),
                                         gos)
            feature['attributes'][attr_key] = attr
        return (FEATURE, feature)
    return go_annot_mapper

def create_dbxref_adder(dbxref_db, relations, dbxref_name='Dbxref'):
    '''It creates a mapper that adds a dbxref to the feature.

    It looks in the provided relations dict to add the dbxref.
    The relations has the feature Name as key and the dbxref as value
    '''
    def add_dbxref_to_feature(item):
        'it adds a dbxref to a feature'
        if item[0] != FEATURE:
            return item
        feature = item[1]
        if feature['id'] in relations:
            dbxref = relations[feature['id']]
            dbxref = dbxref_db + ':' + dbxref
            dbxref = _add_dbxrefs_to_dbxref(feature['attributes'].get(dbxref_name,
                                                                      ''),
                                            [dbxref])
            feature['attributes'][dbxref_name] = dbxref
        return FEATURE, feature
    return add_dbxref_to_feature

def _add_dbxrefs_to_dbxref(dbxref_str, new_dbxrefs):
    '''It adds the dbxrefs to the feature.

    If the dbxref tag is not in feature attributes it creates it'''
    #if the dbxref_db is already in the feature we have to update
    if dbxref_str:
        #'GO:001,GO:002,SO:001'
        dbxrefs = set(dbxref_str.split(','))
    else:
        dbxrefs = set()

    #to the current list of dbxrefs we add one more
    #The final result can be Dbxref:"EMBL:10","EMBL:20"
    #This does not make sense for Dbxref but it does for GO
    dbxrefs = dbxrefs.union(new_dbxrefs)

    return ','.join(dbxrefs)

def modify_gff(ingff3_fpath, outgff3_fpath, mappers=None, filters=None):
    'It modifies the gff features with the given mappers and filters'
    in_gff  = GffFile(fpath=ingff3_fpath)
    out_gff = GffFile(fpath=outgff3_fpath, mode='w')
    if filters is None:
        filters = []
    if mappers is None:
        mappers = []
    for item in in_gff.items:
        for mapper in mappers:
            item = mapper(item)
        result = True
        for filter_ in filters:
            result = filter_(item)
            if not result:
                break
        if not result:
            continue
        if item is not None:
            out_gff.write(item)
    out_gff.flush()

class SeqGffWriter(object):
    'It writes sequences in an gff style'
    def __init__(self, fhand, default_type=None, source='.'):
        'It inits the class'
        self.num_features = 0

        if default_type is None:
            default_type = 'region' # SO:000001
        self._default_type = default_type
        self._source = source
        self._writer = GffFile(fhand.name, mode='w')

    def write(self, sequence):
        'It does the real write of the features'
        seq_feature = self._get_seq_feature(sequence)
        self._writer.write((FEATURE, seq_feature))

        #subfeature
        for feature in self._get_sub_features(sequence):
            self._writer.write((FEATURE, feature))
            self.num_features += 1

    def _get_seq_sofa_type(self, sequence):
        'It gets the type of the feature'
        if 'SO' in sequence.annotations:
            return sequence.annotations['SO']
        else:
            return self._default_type

    def _get_sequence_attributes(self, sequence):
        '''It writes gff attributes looking in features and annotations of the
        sequnce'''
        attributes = {}
        attributes['name'] = sequence.name

        if (sequence.description and
            sequence.description != "<unknown description>"):
            attributes['description'] = sequence.description

        if 'GOs' in sequence.annotations:
            gos = ','.join(sequence.annotations['GOs'])
            attributes['Ontology_term'] = gos

        #orthologs
        orthologs = []
        for annot in sequence.annotations:
            if 'orthologs' in annot:
                specie = annot.split('-')[0]
                for ortholog_names in sequence.annotations[annot]:
                    orthologs.append('%s:%s' % (specie, ortholog_names))
        if orthologs:
            attributes['orthologs'] = ",".join(orthologs)

        return attributes

    def _get_seq_feature(self, sequence):
        'It gets the gff section of the sequence. The parent'
        feature = {}
        feature['seqid'] = sequence.id
        feature['source'] = self._source
        feature['type'] = self._get_seq_sofa_type(sequence)
        feature['start'] = 1
        feature['end'] = len(sequence)

        feature['attributes'] = self._get_sequence_attributes(sequence)

        return feature

    def _get_sub_features(self, sequence):
        'It gets the features of the sequence feature'
        for seq_feature in sequence.features:
            feature = {}
            feature['type'] = seq_feature.type
            feature['seqid'] = sequence.id
            feature['start'] = int(str(seq_feature.location.start)) + 1
            feature['end'] = int(str(seq_feature.location.end)) + 1
            feature['source'] = 'franklin'
            attributes = {}
            if feature['type'] == 'microsatellite':
                #'microsatellite' #SO:0000289
                attributes['score'] = str(seq_feature.qualifiers['score'])
            elif feature['type'] == 'intron':
                feature['type'] = 'intron_' #SO:0000188
            elif feature['type'] == 'orf':
                feature['type'] = 'ORF' #SO:0000236
                strand = seq_feature.qualifiers['strand']
                attributes['strand'] = '+' if strand == 'forward' else '-'
            elif feature['type'] == 'snv':
                feature['type'] = 'SNV' #SO:0001483
            attributes['name'] = feature['seqid'] + '_' + feature['type']
            feature['attributes'] = attributes
            yield feature

    @staticmethod
    def _get_subfeature_attributes(id, name, kind, num):
        '''It gets the attribute section of a sequence's  subfeature'''
        return 'ID=%s_%s_%d;name=%s_%s_%d' % (id, kind, num, name, kind, num)
