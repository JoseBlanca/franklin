'''Code to write cmap gff3 files
Created on 27/10/2009
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
from copy import copy

from biolib.gff import write_gff

def _species_pragma(species):
    'It returns the species pragma line'
    str_ = '##cmap_species '
    str_ += 'species_acc=%s;' % species['accession']
    if 'full_name' in species:
        str_ += 'species_full_name=%s;' % species['full_name']
    if 'common_name' in species:
        str_ += 'species_common_name=%s;'% species['common_name']
    if 'display_order' in species:
        str_ += 'display_order=%s;' % str(species['display_order'])
    return str_

def _map_set_pragma(map_set):
    'It returns the map set pragma line'
    str_ = '##cmap_map_set map_set_acc=%s;' % map_set['accession']
    str_ += 'map_set_name=%s;' % map_set['name']
    if 'short_name' in map_set['short_name']:
        str_ += 'map_set_short_name=%s;' % map_set['short_name']
    str_ += 'map_type_acc=%s;' % map_set['type']
    str_ += 'unit_modifier=%.2f;' % map_set['unit_modifier']
    return str_

def _map_pragma(map_, map_set_accession):
    'It returns the map_ pragma line'
    str_ = '##cmap_map map_acc=%s;' % map_['accession']
    str_ += 'map_name=%s;' % map_['name']
    str_ += 'map_start=%d;' % map_['start']
    str_ += 'map_stop=%d;' % map_['end']
    if 'display_order' in map_:
        str_ += 'display_order=%s;' % str(map_['display_order'])
    str_ += 'map_set_acc=%s;' % map_set_accession
    str_ += '\n'
    str_ += '##sequence-region %s %d %d' % (map_['name'], map_['start'],
                                            map_['end'])
    return str_

def _map_features(map_, features):
    'It returns the list of features for this map'
    feats = []
    for feat_loc in map_['feature_locations']:
        feat = copy(features[feat_loc['feature']])
        feat['seqid'] = map_['name']
        feat['id'] = feat['name']
        feat['source'] = 'CMap'
        feat['start'] = feat_loc['start']
        if 'end' in feat_loc:
            feat['end'] = feat_loc['end']
        else:
            feat['end'] = feat_loc['start']
        feats.append(feat)
    return feats

def cmap_to_gff(data, fhand):
    'Given a dict with the cmap data and an output fhand it writes a gff3 file'
    gff = []
    gff.append('##cmap-gff-version 1')

    for mapset in data['map_sets']:
        species_name = mapset['species']
        species = data['species'][species_name]
        gff.append(_species_pragma(species))
        gff.append('###')
        gff.append(_map_set_pragma(mapset))
        for map_ in mapset['maps']:
            #start and end
            start = None
            end = None
            for feat_loc in map_['feature_locations']:
                this_start = feat_loc['start']
                if 'end' in feat_loc:
                    this_end = feat_loc['end']
                else:
                    this_end = feat_loc['start']
                if start is None or start > this_start:
                    start = this_start
                if end is None or end < this_end:
                    end = this_end
            map_['start'] = start
            map_['end'] = end
            gff.append(_map_pragma(map_, mapset['accession']))
            gff.extend(_map_features(map_, data['features']))


    write_gff(gff, fhand)