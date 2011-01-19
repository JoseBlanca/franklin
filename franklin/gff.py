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

ESCAPES = {';': '%3B',
           '=': '%3D',
           '&': '%26',
           ' ': '%20',
           '%': '%25',
           ',': '%2C'}

def _deescape(string):
    'It replaces the escaped sequences with the real characters'

    for character, escape in ESCAPES.items():
        string = string.replace(escape, character)
    return string

def create_feature_id(feature_name, feature_ids):
    'It returns a new unique feature_id based on the name'
    if feature_name not in feature_ids:
        feature_id = feature_name
        feature_ids[feature_name] = 1
    else:
        feature_id = '%s_%i' % (feature_name, feature_ids[feature_name])

    feature_ids[feature_name] += 1
    return feature_id

def features_in_gff(fhand, version):
    '''It parses a gff file and returns an iterator of each line of the gff
    parsed'''

    feature_ids = {}

    for line in fhand:
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            continue
        yield create_feature(line, version, feature_ids)

def create_feature(line, version, feature_ids=None):
    'It creates a feature from a line of a gff'
    feature = {}
    seqid, source, type_, start, end, score, strand, phase, annots = \
                                                        line.split("\t" , 8)
    start = int(start)
    end   = int(end)

    attributes = {}
    if version == 2:
        separator = ' '
    elif version == 3:
        separator = '='

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
        feature_id = create_feature_id(feature['name'], feature_ids)

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

def get_gff_header(fhand):
    'it gets the gff headers of a gff file'
    original_pos = fhand.tell()
    fhand.seek(0)
    headers = []
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        if not line.startswith('##'):
            fhand.seek(original_pos)
            return headers
        headers.append(line)
    return headers

def add_dbxref_to_feature(feature, dbxref_db, dbxref_id):
    '''It adds the dbxref to the feature. If the dbxref tag is not in feature
    attributes it creates it'''
    dbxref = '%s:%s' % (dbxref_db, dbxref_id)
    if 'Dbxref' in feature['attributes']:
        feature['attributes']['Dbxref'] += ',%s' % dbxref
    else:
        feature['attributes']['Dbxref'] = dbxref

class GffWriter(object):
    'It writes GFF files'
    def __init__(self, fhand, header=None):
        'It inits the class and writes the header'
        self._fhand = fhand

        #write the header
        if header is None:
            fhand.write('##gff-version 3\n')
        else:
            for header_line in header:
                fhand.write('%s\n' % header_line)

        self._feature_ids ={}

    def write(self, feature):
        'It writes a feature in to the gff file'
        if isinstance(feature, dict):
            self._fhand.write(self._feature_to_str(feature))
        else:
            self._fhand.write(feature)
            self._fhand.write('\n')

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
            feature_id = create_feature_id(attributes['Name'],
                                            self._feature_ids)

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

        return '\t'.join(feature_fields) + '\n'

def write_gff(features, out_fhand, header=None):
    '''It writes a gff file.

    It requires a list of features or directives.
    A feature is a dict with the following keys: seqid, source, kind, start,
    end, score, strand, phase, id, name, alias, target, gap, derives_from, note,
    dbxred, ontology_term, parents, and children.
    A directive is just a plain str.
    '''
    writer = GffWriter(fhand=out_fhand, header=header)

    for feature in features:
        writer.write(feature)

