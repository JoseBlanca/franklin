'''This module deals with the IO of GFF files.

Created on 26/10/2009
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

def _feature_to_str(feature):
    'Given a feature dict it returns a gff feature line'
    feat_str = []

    feat_str.append(feature['seqid'])
    feat_str.append('\t')

    if 'source' in feature:
        feat_str.append(feature['source'])
    else:
        feat_str.append('.')
    feat_str.append('\t')

    feat_str.append(feature['kind'])
    feat_str.append('\t')

    feat_str.append('%d' % feature['start'])
    feat_str.append('\t')

    feat_str.append('%d' % feature['end'])
    feat_str.append('\t')

    for tag, default in (('score', '.'), ('strand', '+'), ('phase', '.')):
        if tag in feature:
            feat_str.append(feature[tag])
        else:
            feat_str.append(default)
        feat_str.append('\t')

    #the attributes
    attributes = []

    if 'id' in feature:
        attributes.append('ID=%s' % feature['id'])
    if 'name' in feature:
        attributes.append('Name=%s' % feature['name'])
    if 'parents' in feature:
        attributes.append('Parent=%s' % ','.join(feature['parents']))

    #the rest of the attributes
    done_attrs = ('seqid', 'source', 'kind', 'start', 'end', 'id', 'name',
                  'parents')
    for attribute, value in feature.items():
        if attribute not in done_attrs:
            attributes.append('%s=%s' % (attribute, value))

    feat_str.append(';'.join(attributes))

    feat_str.append('\n')
    return ''.join(feat_str)

def write_gff(features, out_fhand):
    '''It writes a gff file.

    It requires a list of features or directives.
    A feature is a dict with the following keys: seqid, source, kind, start,
    end, score, strand, phase, id, name, alias, target, gap, derives_from, note,
    dbxred, ontology_term, parents, and children.
    A directive is just a plain str.
    '''
    out_fhand.write('##gff-version 3\n')

    for feat in features:
        if isinstance(feat, dict):
            out_fhand.write(_feature_to_str(feat))
        else:
            out_fhand.write(feat)
            out_fhand.write('\n')
