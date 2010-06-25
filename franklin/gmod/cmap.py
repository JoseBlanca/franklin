'''Code to write cmap gff3 files
Created on 27/10/2009
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
from copy import copy

from itertools import combinations
from franklin.gff import write_gff

def _species_pragma(species):
    'It returns the species pragma line'
    str_ = '##cmap_species '
    str_ += 'species_acc=%s;' % species['accession']
    if 'full_name' in species:
        str_ += 'species_full_name=%s;' % species['full_name']
    if 'common_name' in species:
        str_ += 'species_common_name=%s;' % species['common_name']
    if 'display_order' in species:
        str_ += 'display_order=%s;' % str(species['display_order'])
    return str_

def _map_set_pragma(map_set):
    'It returns the map set pragma line'
    str_ = '##cmap_map_set map_set_acc=%s;' % map_set['accession']
    str_ += 'map_set_name=%s;' % map_set['name']
    if 'short_name' in map_set['short_name']:
        str_ += 'map_set_short_name=%s;' % map_set['short_name']
    else:
        str_ += 'map_set_short_name=%s;' % map_set['accession']
    str_ += 'map_type_acc=%s;' % map_set['type']
    str_ += 'unit_modifier=%.2f;' % map_set['unit_modifier']
    return str_

def _map_pragma(map_, map_set_accession):
    'It returns the map_ pragma line'
    str_ = '##cmap_map map_acc=%s_%s;' % (map_set_accession, map_['accession'])
    str_ += 'map_name=%s_%s;' % (map_set_accession, map_['name'])
    str_ += 'map_start=%d;' % (map_['start'])
    str_ += 'map_stop=%d;' % map_['end']
    if 'display_order' in map_:
        str_ += 'display_order=%s;' % str(map_['display_order'])
    str_ += 'map_set_acc=%s;' % map_set_accession
    str_ += '\n'
    str_ += '##sequence-region %s_%s %d %d' % (map_set_accession, map_['name'],
                                               map_['start'], map_['end'])
    return str_

def _map_features(map_, features, map_set_accession, marker_count,
                  marker_id_map, features_in_mapset):
    'It returns the list of features for this map'
    feats = []

    for feat_loc in map_['feature_locations']:
        feat = copy(features[feat_loc['feature']])
        feat['seqid'] = '%s_%s' % (map_set_accession, map_['name'])
        marker_name = feat['name']
        feat_acc = (feat['seqid'], marker_name)
        if  feat_acc not in features_in_mapset:
            features_in_mapset.add(feat_acc)
        else:
            raise RuntimeError('%s marker already in %s mapset' % \
                                            (marker_name, feat['seqid']))

        if marker_name not in marker_count:
            marker_count[marker_name] = []

        marker_id = '%s_%d' % (marker_name, len(marker_count[marker_name]))
        marker_count[marker_name].append(marker_id)

        # construction of marker id map
        marker_id_map[marker_id] = (map_set_accession, map_['accession'])

        feat['id'] = marker_id
        feat['source'] = 'CMap'
        feat['start'] = feat_loc['start']
        if 'end' in feat_loc:
            feat['end'] = feat_loc['end']
        else:
            feat['end'] = feat_loc['start']
        if 'attributes' not in feat:
            feat['attributes'] = {}
        feat['attributes']['feature_accs'] = '%s_%s' %  (feat['seqid'],
                                                         marker_name)

        feats.append(feat)
    return feats

def _cmap_correspondences(marker_count, marker_id_map):
    '''It prints the correspondence pragmas. It needs a list of all the
    correspondences'''
    strs = []
    for correspondences in marker_count.values():
        for marker1, marker2 in combinations(correspondences, 2):
            str_ = ['##cmap_corr evidence_type_acc=ANB;']
            str_.append('ID1=%s;map_set_acc1=%s;map_acc1=%s;' % (marker1,
                                                 marker_id_map[marker1][0],
                                                 marker_id_map[marker1][1]))
            str_.append('ID2=%s;map_set_acc2=%s;map_acc2=%s;' % (marker2,
                                                 marker_id_map[marker2][0],
                                                 marker_id_map[marker2][1]))
            strs.append(''.join(str_))
    return strs

def cmap_to_gff(data, fhand):
    'Given a dict with the cmap data and an output fhand it writes a gff3 file'
    gff = []
    gff.append('##cmap-gff-version 1')
    # This marker count is used where there is a  markers in two maps.
    marker_count = {}
    marker_id_map = {}
    for mapset in data['map_sets']:
        species_name = mapset['species']
        species = data['species'][species_name]
        gff.append(_species_pragma(species))
        gff.append('###')
        gff.append(_map_set_pragma(mapset))
        features_in_mapset = set()
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
            gff.extend(_map_features(map_, data['features'],
                                     mapset['accession'], marker_count,
                                     marker_id_map,
                                     features_in_mapset))
    #the correspondences
    gff.extend(_cmap_correspondences(marker_count, marker_id_map))

    write_gff(gff, fhand)

def cmap_to_mcf(data, fhand):
    'Given a dict with the cmap data and an output fhand it writes a mcf file'

    for mapset in data['map_sets']:
        fhand.write('**************************************************\n')
        fhand.write('map: %s\n' % mapset['name'])
        fhand.write('**************************************************\n\n')
        for map_ in mapset['maps']:
            fhand.write(map_['name'] + '\n')
            markers = []
            qtls = []
            for feat_loc in map_['feature_locations']:
                start = feat_loc['start']
                if 'end' in feat_loc:
                    end = feat_loc['end']
                else:
                    end = feat_loc['start']
                if start == end:
                    markers.append({'name':feat_loc['feature'],
                                    'start':start})
                else:
                    qtls.append({'name':feat_loc['feature'],
                                 'start':start,
                                 'end':end})
            #now we sort the markers
            numeric_sort = lambda x, y: int(x['start'] - y['start'])
            markers = sorted(markers, numeric_sort)
            for marker in markers:
                fhand.write('%s\t%s\n' % (marker['name'], marker['start']))
            fhand.write('qtls\n')
            qtls = sorted(qtls, numeric_sort)
            for marker in qtls:
                location = float(marker['start'] + marker['end']) / 2.0
                location = '%.1f' % location
                fhand.write('%s\t%s\t%s\t%s\t%s\n' % (marker['name'], location,
                                                    location, location,
                                                    location))
        fhand.write('\n\n')
