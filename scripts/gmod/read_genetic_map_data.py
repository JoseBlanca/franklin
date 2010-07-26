#!/usr/bin/python
'It reads the genetic map data from the csv file'

#pylint: disable-msg=C0301

import csv, os.path
from franklin.gmod.cmap import cmap_to_gff, cmap_to_mcf
from franklin.gmod.markers import (parse_markersfile,
                                   get_marker_from_correlations)

def get_correlated_name(name, correlations):
    'It return the name fixed by the correlation'
    # it looks in the marker_correlations_file
    marker_ =  get_marker_from_correlations(name, correlations)
    if marker_ is not None:
        name = marker_['name']
    return name

def read_icugi_map(fhand, cmap, correlations=None):
    'It reads the merged icugi map from the text file'
    icugi = 'icugi'
    #the mapset
    cmap['map_sets'].append({'species':'cmelo', 'accession':icugi,
                             'name':icugi,'short_name':icugi,
                             'type':'Genetic', 'unit_modifier': 0.1,
                             'maps':[]})
    #map_set_index
    mapset_index = {}
    for index, mapset in enumerate(cmap['map_sets']):
        mapset_index[mapset['name']] = index

    map_index = {}

    markers = cmap['features']
    icugi_mapset = cmap['map_sets'][mapset_index[icugi]]
    group = None
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        #the group
        if line.startswith('group'):
            group = line.split()[1].split('_')[1].lower()
            #print group
            icugi_mapset['maps'].append({'accession': group,
                                          'name': group,
                                          'display_order':1,
                                          'feature_locations':[]})
            map_index[group] = len(icugi_mapset['maps']) - 1
            continue
        line_items = line.split()
        marker, location = line_items[0], line_items[1]
        location = location.replace(',', '.')

        location = int(round(float(location)))

        #the marker
        marker = marker.lower()
        marker = marker.strip()
        # it looks in the marker_correlations_file
        marker = get_correlated_name(marker, correlations)

        if marker not in markers:
            #print markers.keys()
            markers[marker] = {'name':marker, 'type':'unknown'}
            print 'marker not found ->', marker

        #the marker in the map
        feat_loc = {'feature':marker, 'start':location}
        group_index = map_index[group]
        icugi_mapset['maps'][group_index]['feature_locations'].append(feat_loc)

def read_maps(fhand, cmap, correlations=None):
    'Given a csv it reads the map and marker data'
    map_reader = csv.reader(fhand, delimiter=',')

    # unit modifier
    unit_modifier = 0.1
    #a index for the columns tags
    cols = {}
    for index, col in enumerate(map_reader.next()):
        cols[col] = index

    map_tree = {'oliver'   :{'maps':{}},
                'gonzalo'  :{'maps':{}},
                'fernandez':{'maps':{}},
                'deleu'    :{'maps':{}},
               }

    markers = cmap['features']
    for row in map_reader:
        group = row[cols['GL new']]
        if not group:
            continue
        name = row[cols['Marker']].lower()
        name = name.strip()
        name = get_correlated_name(name, correlations)

        if row[cols['Type']] == 'SNP-indel, RFLP':
            type_ = 'SNP-INDEL'
        else:
            type_ = row[cols['Type']]
        marker = {'name':name,
                  'type':type_,
                  'publication':row[cols['Published']],
                 }
        alias = row[cols['Other names']]
        if alias:
            marker['alias'] = alias
        est = row[cols['EST']]
        if est:
            marker['est'] = est

        markers[name] = marker

        # creation and population of the maps an the feature locs

        for map_, col in (('deleu', 'Deleu Bin'), ('oliver', 'Oliver F2'),
                          ('gonzalo', 'Gonzalo DHL'),
                          ('fernandez', 'Fernandez DHL'),
                          ('unpub', 'Unpublished Bin')):
            position = row[cols[col]]
            if position:
                position = float(position.replace(',', '.'))
                position = position / float(unit_modifier)
                if map_ == 'unpub':
                    map_ = 'deleu'
                if group not in map_tree[map_]['maps']:
                    map_tree[map_]['maps'][group] = {'accession': group.lower(),
                                                     'name': group.lower(),
                                                     'display_order':1,
                                                     'feature_locations':[]}
                feat_loc = {'feature':name, 'start':position}
                map_tree[map_]['maps'][group]['feature_locations'].append(feat_loc)

    #Create maps set structure
    mapsets = {'oliver':{'species':'cmelo', 'accession':'oliver',
                         'name':'oliver','short_name':'oliver',
                         'type':'Genetic', 'unit_modifier': unit_modifier,
                         'maps':map_tree['oliver']['maps'].values()},
               'gonzalo':{'species':'cmelo', 'accession':'gonzalo',
                          'name':'gonzalo', 'short_name':'gonzalo',
                          'type':'Genetic', 'unit_modifier': unit_modifier,
                          'maps':map_tree['gonzalo']['maps'].values()},
               'fernandez': {'species':'cmelo', 'accession':'fernandez',
                          'name':'fernandez', 'short_name':'fernandez',
                          'type':'Genetic',  'unit_modifier': unit_modifier,
                          'maps': map_tree['fernandez']['maps'].values()},
               'deleu': {'species':'cmelo', 'accession':'deleu', 'name':'deleu',
                         'short_name':'deleu', 'type':'Genetic',
                         'unit_modifier': unit_modifier,
                         'maps':map_tree['fernandez']['maps'].values()}}
    cmap['map_sets'].extend(mapsets.values())
    #pprint.pprint(cmap)

def create_empty_cmap():
    'It returns an empty cmap structure'
    # Species is hard coded.
    species = {'accession':'cmelo',
               'full_name':'Cucumis melo',
               'common_name':'melon',
               'display_order':1}
    cmap = {'features':{},
            'map_sets':[],
            'species':{'cmelo': species}}
    return cmap

def read_markers(fhand, cmap, marker_corr=None):
    '''It reads the Monforte marker data. It uses the correlations file to
    adjust data'''
    map_reader = csv.reader(fhand, delimiter='\t')

    markers = cmap['features']
    for row in map_reader:
        name = row[0].lower()
        name = name.strip()
        kind = row[1].strip()
        reference = row[2]
        correlation = get_marker_from_correlations(name, marker_corr)
        if correlation is not None:
            name = correlation['name']
            if correlation['type']  is None:
                kind = correlation['sofa']
            else:
                kind = correlation['type']
            reference = correlation['publication']
        if not kind:
            kind = 'unknown'
        markers[name] = {'name':name,
                         'type':kind,
                          'publication':reference}

def fix_marker_types(cmap):
    'It standarizes the marker types'

    std_types = {
        'aflp': ('aflp', ),
        'est-ssr' : ('microsatellite', 'SO:0000289'),
        'ssr': ('microsatellite', 'SO:0000289'),
        'ima': ('issr', ),
        'indel': ('indel', 'SO:1000032'),
        'snp-indel': ('indel', 'SO:1000032'),
        'marker' : ('sequence_feature', 'SO:0000110'),
        'morhological': ('morphological', ),
        'morphological': ('morphological', ),
        'trait': ('morphological', ''),
        'rapd': ('rapd', ),
        'rapds': ('rapd',),
        'rflp': ('RFLP_fragment', ),
        'snp': ('SNP', 'SO:0000694'),
        'snp-caps': ('SNP', 'SO:0000694', 'caps'),
        'snp-snapshot-caps': ('SNP', 'SO:0000694', 'caps'),
        'snp-scar': ('SNP', 'SO:0000694', 'caps'),
        'scar': ('SNP', 'SO:0000694', 'caps'),
        'snp-sequencing': ('SNP', 'SO:0000694'),
        'snp-snapshot': ('SNP', 'SO:0000694'),
        'spelling error': ('spelling error',),
        'unknown': ('sequence_feature', 'SO:0000110'),
        'isozyme' : ('isozyme', )
    }

    markers = cmap['features']
    for marker in markers:
        type_ = markers[marker]['type'].lower()
        if type_ in std_types:
            type_ = std_types[type_][0]
            markers[marker]['type'] = type_

class Cmap(object):
    'Some useful function for a map set'
    def __init__(self, cmap):
        'it requires a cmap structure'
        self.cmap = cmap
        self._map_index = {}
        self._index = {}

    def _create_indexes(self):
        'It creates an index for a mapset'
        #first we create the mapset index
        index = {}
        for mapset_index, mapset in enumerate(self.cmap['map_sets']):
            mapset = mapset['accession']
            for map_ in self.cmap['map_sets'][mapset_index]['maps']:
                group = map_['accession']
                #print mapset, group
                if mapset not in index:
                    index[mapset] = {}
                index[mapset][group] = map_
        self._map_index = index
        #now we create the marker position index
        for mapset in self._map_index:
            for map_name, map_ in self._map_index[mapset].items():
                for feature_location in map_['feature_locations']:
                    #print 'mapset', mapset, map_name, feature_location
                    if mapset not in self._index:
                        self._index[mapset] = {}
                    if map_name not in self._index[mapset]:
                        self._index[mapset][map_name] = {}
                    #print mapset, map_name, feature_location['feature']
                    self._index[mapset][map_name][feature_location['feature']] = \
                                                            feature_location
    def feature_location(self, feature, mapset, map_=None):
        'It returns the position of  a marker in a map'
        if mapset not in self._map_index:
            self._create_indexes()
        feat_loc = {'start':None, 'end':None}
        if map_:
            try:
                #print self._index[mapset][map_][feature]
                feat_loc = self._index[mapset][map_][feature]
            except KeyError:
                raise KeyError('%s feature not in %s map' % (feature, map_))

        else:
            for map_ in self._index[mapset]:

                if feature in self._index[mapset][map_]:
                    feat_loc = self._index[mapset][map_][feature]
                    break
            else:
                #print feature, self._index[mapset]
                raise ValueError('Marker not found mapset %s feature %s' % (mapset, feature))


        start = feat_loc['start']
        if 'end' in feat_loc:
            end = feat_loc['end']
        else:
            end = start
        return mapset, map_, start, end

    def get_map(self, mapset, map_):
        'It returns a map'
        if mapset not in self._map_index:
            self._create_indexes()
        return self._map_index[mapset][map_]

def read_bin_map(fhand, cmap, correlations):
    'It loads the bin map '
    map_reader = csv.reader(fhand, delimiter='\t')
    map_reader.next() #ignore fist line

    map_data = []
    for row in map_reader:
        group       = row[1]
        marker      = row[2].lower().strip()
        marker_type = row[3]
        bin_        = row[5]
        if not group and not bin_:
            continue
        marker_name = get_correlated_name(marker, correlations)
        marker_ =  get_marker_from_correlations(marker_name, correlations)
        feature_type = marker_['sofa']

        data = { 'group': group.lower(), 'marker':marker_name,
                        'marker_type': marker_type, 'bin': bin_,
                        'feature_type':feature_type}
        map_data.append(data)


    all_bins  = []
    list_bins = []
    ref       = None
    group     = None
    for marker in map_data:
        kind        =  marker['marker_type']
        marker_name = marker['marker']

        if group is None or group !=  marker['group']:
            for bin_marker in list_bins:
                bin_marker['end_marker']   = ref
                all_bins.append(bin_marker)
                list_bins = []
            ref = None
        group = marker['group']

        if kind in ('bin', 'binint'):
            bin_marker = {'start_marker': None, 'end_marker': None,
                          'type':kind, 'name':marker_name, 'map': group,
                          'feature_type':marker['feature_type']}
            if ref is None:
                bin_marker['start_marker'] = None
            else:
                bin_marker['start_marker'] = ref
            list_bins.append(bin_marker)

        elif kind == 'ref':
            for bin_marker in list_bins:
                bin_marker['end_marker']   = marker_name
                all_bins.append(bin_marker)
            list_bins = []
            ref = marker_name

    cmap_ = Cmap(cmap)
    # now put the positions
    for bin in all_bins:
        if bin['start_marker'] is None:
            bin['start_position'] = 0
        else:
            try:
                start = cmap_.feature_location(feature=bin['start_marker'],
                                                    mapset='icugi',
                                                    map_=bin['map'])[2]

            except KeyError:
                print "bin start not found, reference: %s group : %s " % \
                                                    (bin['start_marker'], bin['map'])
                start = 0 # this is to be able to work, this is unreal
            bin['start_position'] = start

        if bin['end_marker'] is None:
            bin['end_postion'] = bin['start_position']
        else:
            try:
                start = cmap_.feature_location(feature=bin['end_marker'],
                                                    mapset='icugi',
                                                    map_=bin['map'])[2]
            except KeyError:
                print "bin end not found, reference: %s group : %s " % \
                                                    (bin['end_marker'], bin['map'])
                start = 0 # this is to be able to work, this is unreal
            bin['end_position'] = start

    # now look at the start position of the bininit
    for bin in all_bins:
        if bin['type'] == 'binint':
            half = (bin['end_position'] - bin['start_position']) / float(2)
            bin['start_position'] += half
            bin['start_position'] = int(round(bin['start_position']))

    #add to cmap
    for bin in all_bins:
        #add to markers
        marker_name = bin['name']
        sofa = bin['feature_type']
        cmap['features'][marker_name] = {'type':sofa, 'name':marker_name}
        # add to mapset(icugi)
        feat_loc = {'feature':marker_name, 'start':bin['start_position'],
                    'end':bin['end_position']}
        # Now we need to guess where is this feature_loc adding.
        # first traduce from icugi to map number
        for index, mapset in enumerate(cmap['map_sets']):
            if mapset['name'] == 'icugi':
                icugi_index = index
                for index2, map_ in enumerate(mapset['maps']):
                    if map_['accession'] == bin['map']:
                        map_index = index2
        cmap['map_sets'][icugi_index]['maps'][map_index]['feature_locations'].append(feat_loc)

def read_qtls(fhand, cmap, correlations=None):
    'It reads the qtls markers for the icugi map'
    def _int(an_str):
        'It converts to int taking into account Nones'
        if not an_str:
            return None
        else:
            return int(an_str)
    qtls = {}
    rows =  csv.reader(fhand, delimiter='\t')
    rows.next() # we don't need first line
    for row in rows:
        qtl_id = row[3].strip().lower()
        if not qtl_id:
            continue
        start_marker = get_correlated_name(row[6].lower().strip(), correlations)
        end_marker   = get_correlated_name(row[7].lower().strip(), correlations)
        qtls[qtl_id] = {'start_marker': start_marker,
                        'end_marker': end_marker,
                        'start_position': None,
                        'end_position': None}
    cmap_ = Cmap(cmap)
    #now we need the position for every marker in the qtls
    for qtl_id, qtl in qtls.items():
        start_marker, end_marker = qtl['start_marker'], qtl['end_marker']
        start_map, end_map = None, None
        start_position, end_position = None, None
        if start_marker:
            try:
                result         = cmap_.feature_location(feature=start_marker,
                                                        mapset='icugi')
                start_map      = result[1]
                start_position = result[2]
            except ValueError:
                print 'for qtl %s start marker not found %s' % (qtl_id,
                                                                start_marker)
        if end_marker:
            try:
                result       = cmap_.feature_location(feature=end_marker,
                                                      mapset='icugi')
                end_map      = result[1]
                end_position = result[3]
            except ValueError:
                print 'for qtl %s end marker not found %s' % (qtl_id,
                                                              end_marker)

        if start_map and end_map:
            #assert start_map == end_map
            if start_map != end_map:
                print '%s %s not in the same map %s %s' % (start_marker,
                                                          end_marker, start_map,
                                                          end_map)

        map_ = start_map
        if start_position is not None:
            qtl['start_position'] = start_position

        if end_position is not None:
            qtl['end_position'] = end_position

        qtl['map']   = map_
        qtls[qtl_id] = qtl


    for qtl_id, qtl_data in qtls.items():
        #add to markers
        cmap['features'][qtl_id] = {'type':'QTL', 'name':qtl_id}
        # add to mapset(icugi)
        feat_loc = {'feature':qtl_id, 'start':qtl_data['start_position'],
                    'end':qtl_data['end_position']}
        # Now we need to guess where is this feature_loc adding.
        # first traduce from icugi to map number
        for index, mapset in enumerate(cmap['map_sets']):
            if mapset['name'] == 'icugi':
                icugi_index = index
                for index2, map_ in enumerate(mapset['maps']):
                    if map_['accession'] == qtl_data['map']:
                        map_index = index2
        cmap['map_sets'][icugi_index]['maps'][map_index]['feature_locations'].append(feat_loc)

def main():
    'It reads the map data and writes it'
    base_path       = '/home/peio/melonomics/datos_originales/mapa_genetico'
    markers_path    = 'mapa_icugi/melon_markers.csv'
    merged_map_path = 'mapa_icugi/framework_merged_map.csv'
    irta_maps_path  = 'mapa_genetico_irta.csv'
    icugi_qtls      = 'mapa_icugi/all_qtls.csv'

    bin_map_path    = 'mapa_icugi/bin_mapping.csv'

    out_path        = 'out.gff3'
    mcf_path        = 'out.mcf'


    correl_fpath    = os.path.join(base_path, 'mapa_icugi/markers.dat')
    markers_path    = os.path.join(base_path, markers_path)
    merged_map_path = os.path.join(base_path, merged_map_path)
    irta_maps_path  = os.path.join(base_path, irta_maps_path)
    out_path        = os.path.join(base_path, out_path)
    mcf_path        = os.path.join(base_path, mcf_path)
    bin_map_path    = os.path.join(base_path, bin_map_path)
    qtls_path       = os.path.join(base_path, icugi_qtls)
    #qtls_fix_path   = os.path.join(base_path, icugi_fix_qtls)

    cmap = create_empty_cmap()

    #Correlations file
    marker_corr = parse_markersfile(open(correl_fpath))

    #get Monforte marker data
    read_markers(open(markers_path), cmap, marker_corr)

    #get irta maps
    cmap_fhand = open(irta_maps_path)
    read_maps(cmap_fhand, cmap, marker_corr)

    #get data from icugi map
    read_icugi_map(open(merged_map_path), cmap, marker_corr)

    #fix the marker types
    fix_marker_types(cmap)


    #get icugi qtls
    read_qtls(open(qtls_path), cmap, marker_corr)

    #read bin map
    read_bin_map(open(bin_map_path), cmap, marker_corr)
#    import pprint
#    pprint.pprint(cmap)
#    return

    #write gff3 file
    gff_fhand  = open(out_path, 'w')
    cmap_to_gff(cmap, gff_fhand)

    #write mcf file
    mcf_fhand  = open(mcf_path, 'w')
    cmap_to_mcf(cmap, mcf_fhand)

if __name__ == '__main__':
    main()

