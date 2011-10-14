'''
Created on 14/06/2010

@author: jose
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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
from __future__ import division
from franklin.seq.seqs import UNKNOWN_DESCRIPTION
from franklin.snv.snv_annotation import calculate_snv_kind, SNV_TYPES
from franklin.statistics import draw_stacked_columns

def _location_to_orf(orfs, feat):
    'It returns the location of the feature respect the orf'
    orfs = list(orfs)
    if not orfs or len(orfs) > 1:
        return 'unknown'

    _loc_to_int = lambda x: (int(x.start.position), int(x.end.position))

    orf_loc = _loc_to_int(orfs[0].location)
    feat_loc = _loc_to_int(feat.location)

    #do they overlap?
    if (feat_loc[0] >= orf_loc[0] and feat_loc[0] <= orf_loc[1] or
        feat_loc[1] >= orf_loc[0] and feat_loc[1] <= orf_loc[1]):
        return 'orf'
    #is feat before orf?
    in_5_prime = True if feat_loc[0] < orf_loc[0] else False
    #is the orf reversed?
    if orfs[0].qualifiers['strand'] != 'forward':
        in_5_prime = not (in_5_prime)
    if in_5_prime:
        return 'utr5'
    else:
        return 'utr3'

MICROSATELLITE_TYPES = ['dinucleotide',
                        'trinucleotide',
                        'tetranucleotide',
                        'pentanucleotide',
                        'hexanucleotide']

def _do_ssr_stats(stats, feats, orfs):
    'It adds the ssr stats'

    types = dict([(index+2, type_) for index, type_ in enumerate(MICROSATELLITE_TYPES)])

    some_feat = False
    for feat in feats:
        some_feat = True
        unit = feat.qualifiers['unit']
        if unit not in stats['microsatellites']['units']:
            stats['microsatellites']['units'][unit] = 0
        stats['microsatellites']['units'][unit] += 1

        type_ = types[len(unit)]
        if type_ not in stats['microsatellites']['types']:
            stats['microsatellites']['types'][type_] = {}
            stats['microsatellites']['types'][type_]['num'] = 0
            stats['microsatellites']['types'][type_]['location'] = {}

        stats['microsatellites']['types'][type_]['num'] += 1
        location = _location_to_orf(orfs, feat)

        if location:
            if location not in stats['microsatellites']['types'][type_]['location']:
                stats['microsatellites']['types'][type_]['location'][location] = 0
            stats['microsatellites']['types'][type_]['location'][location] += 1

    if some_feat:
        stats['microsatellites']['n_seqs'] += 1

def _do_snv_stats(stats, feats, orfs):
    'It adds the ssr stats'

    some_feat = False
    for feat in feats:
        some_feat = True
        type_ = calculate_snv_kind(feat, detailed=True)
        if type_ not in stats['snvs']['types']:
            stats['snvs']['types'][type_] = 0
        stats['snvs']['types'][type_] += 1

        location = _location_to_orf(orfs, feat)
        if location:
            if location not in stats['snvs']['locations']:
                stats['snvs']['locations'][location] = 0
            stats['snvs']['locations'][location] += 1
        stats['snvs']['n_snvs'] += 1
    if some_feat:
        stats['snvs']['n_seqs'] += 1


def do_annotation_statistics(seqs, out_fhand):
    'Given some seqs it calculates the annotation statistics'
    stats = _calculate_annot_stats(seqs)
    _write_annot_stats(stats, out_fhand)

def _write_snp_annot_stats(stats, out_fhand):
    'It writes the snp annot stats to a file'
    stats = stats['snvs']
    if not stats['n_seqs']:
        return

    out_fhand.write('\nSNVs\n')
    out_fhand.write('____\n')

    out_fhand.write('Sequences with SNVs: %i\n' % stats['n_seqs'])
    out_fhand.write('SNVs found: %i\n' % stats['n_snvs'])
    out_fhand.write('SNV types:\n')
    for type_ in sorted(stats['types']):
        str_type = SNV_TYPES[type_]
        out_fhand.write('\t%s: %i\n' % (str_type, stats['types'][type_]))

    out_fhand.write('SNV locations:\n')
    for type_ in ('in 5 prime', 'in orf', 'in 3 prime', 'unknown'):
        if type_ in stats['locations']:
            out_fhand.write('\t%s: %i\n' % (type_, stats['locations'][type_]))

def _write_ssr_annot_stats(stats, out_fhand):
    'It writes the ssr annot stats to a file'
    stats = stats['microsatellites']

    if not stats['n_seqs']:
        return

    out_fhand.write('\nMicrosatellites\n')
    out_fhand.write('_______________\n')

    out_fhand.write('Sequences with microsatellites: %i\n' % stats['n_seqs'])

    header = '%8s|%15s|%15s|%15s|%15s|%15s|%8s|' % (' ',
                                                      MICROSATELLITE_TYPES[0],
                                                      MICROSATELLITE_TYPES[1],
                                                      MICROSATELLITE_TYPES[2],
                                                      MICROSATELLITE_TYPES[3],
                                                      MICROSATELLITE_TYPES[4],
                                                      'Total')
    out_fhand.write(header + '\n')
    out_fhand.write('-' * len(header) + '\n' )


    for loc in ('utr3', 'utr5', 'orf', 'unknown'):
        out_fhand.write('%8s|' % loc)
        loc_sum = 0
        for type_ in MICROSATELLITE_TYPES:
            try:
                ssr_stat = stats['types'][type_]['location']

                value = ssr_stat.get(loc, 0)
            except KeyError:
                value = 0

            out_fhand.write('%15s|' % (value))
            loc_sum += value
        out_fhand.write('%8s|\n' % loc_sum)

    out_fhand.write('-' * len(header) + '\n' )
    out_fhand.write('%8s|' % 'total')

    for type_ in MICROSATELLITE_TYPES:
        try:
            value = stats['types'][type_]['num']
        except KeyError:
            value = 0
        out_fhand.write('%15s|' % value)
    out_fhand.write('%8s|\n' % ' ')




def _write_annot_stats(stats, out_fhand):
    'It writes the annot stats to a file'

    msg = 'Annotation statistics'
    out_fhand.write(msg + '\n')
    out_fhand.write('-' * len(msg) + '\n')

    out_fhand.write('Number of sequences: %i\n' % stats['total_seqs'])
    out_fhand.write('Sequences with description: %i\n' %
                                                        stats['seqs_with_desc'])

    if stats['orf']['n_seqs'] != 0:
        out_fhand.write('Sequences with ORF: %i\n' % stats['orf']['n_seqs'])
        out_fhand.write('Number of ORFs: %i\n' % stats['orf']['n_feats'])
    if stats['intron']['n_seqs'] != 0:
        out_fhand.write('Sequences with intron: %i\n' %
                                                     stats['intron']['n_seqs'])
        out_fhand.write('Number of introns: %i\n' % stats['intron']['n_feats'])
    if stats['orthologs']:
        out_fhand.write('\nOrthologs\n')
        out_fhand.write('_________\n')
        for db_ in stats['orthologs']:
            out_fhand.write('Sequences with %s orthologs: %i\n' %
                                       (db_, stats['orthologs'][db_]['n_seqs']))
            out_fhand.write('Number of %s orthologs: %i\n' %
                                     (db_, stats['orthologs'][db_]['n_annots']))
    if stats['GOs']['n_seqs']:
        out_fhand.write('\nGO terms\n')
        out_fhand.write('________\n')

        out_fhand.write('Sequences with GOs: %i\n' % stats['GOs']['n_seqs'])
        out_fhand.write('Number of GOs: %i\n' % stats['GOs']['n_annots'])

    _write_snp_annot_stats(stats, out_fhand)
    _write_ssr_annot_stats(stats, out_fhand)

    out_fhand.write('\n')

def _calculate_annot_stats(seqs):
    'Given some seqs it calculates the annotation statistics'

    annot_stats = {'total_seqs': 0,
                   'seqs_with_desc':0,
                   'orf': {'n_seqs':0,
                           'n_feats':0,},
                   'intron': {'n_seqs':0,
                              'n_feats':0,},
                   'microsatellites': {'n_seqs':0,
                                       'units':{},
                                       'types':{},
                                       'locations':{},},
                   'snvs': {'n_snvs':0,
                            'n_seqs':0,
                            'types':{},
                            'locations':{},},
                   'orthologs':{},
                   'GOs':{'n_seqs':0, 'n_annots':0},
                   }
    for seq in seqs:
        #description
        desc = seq.description
        if desc and desc != UNKNOWN_DESCRIPTION:
            annot_stats['seqs_with_desc'] += 1

        #features
        for feat_kind in ('orf', 'intron'):
            feat_found = False
            for feat in seq.get_features(feat_kind):
                annot_stats[feat_kind]['n_feats'] += 1
                feat_found = True
            if feat_found:
                annot_stats[feat_kind]['n_seqs'] += 1

        #microsatellites and snvs
        orfs = seq.get_features('orf')
        _do_ssr_stats(annot_stats, seq.get_features('microsatellite'),
                      orfs)
        _do_snv_stats(annot_stats, seq.get_features('snv'),
                      orfs)

        #annotations
        for annot, value in seq.annotations.items():
            if 'ortholog' in annot:
                annot_type = 'orthologs'
                annot = annot.split('-')[0]
                if annot not in annot_stats['orthologs']:
                    annot_stats[annot_type][annot] = {'n_seqs':0, 'n_annots':0}
                annot_stats[annot_type][annot]['n_annots'] += len(value)
                annot_stats[annot_type][annot]['n_seqs'] += 1
            elif 'GOs' == annot:
                annot_stats[annot]['n_annots'] += len(value)
                annot_stats[annot]['n_seqs'] += 1
            else:
                continue

        annot_stats['total_seqs'] += 1
    return annot_stats

def _nucleotide_freq_per_position(sequences, positions_to_study=30):
    'It calculates the probability of each nucleotide in each'
    pos_count = []

    for index in range(positions_to_study):
        pos_count.append({'A':0, 'T':0, 'C':0, 'G':0})

    for sequence in sequences:
        positions = positions_to_study
        sequence_seq = str(sequence.seq)
        if len(sequence_seq) < positions_to_study:
            positions = len(sequence_seq)
        for index in range(positions):
            nucl = sequence_seq[index]
            nucl = nucl.upper()
            if nucl not in ('A', 'T', 'C', 'G'):
                continue
            pos_count[index][nucl] += 1

    freq_stats = {'A':[], 'T':[], 'C':[], 'G':[]}
    for index, position in enumerate(pos_count):
        total_nucl = sum(position.values())
        if total_nucl < 0.001:  #from this position onwards there's no stats
            break
        for nucl, count in position.items():
            if total_nucl == 0:
                freq = 0
            else:
                freq = count/total_nucl

            freq_stats[nucl].append(freq)
    return freq_stats

def create_nucleotide_freq_histogram(sequences, fhand=None, title=None,
                                     positions_to_study=30):
    '''It writes a especific stacked_columns graphic representing the freq of
    each nucleotide per position'''
    values = _nucleotide_freq_per_position(sequences,
                                          positions_to_study=positions_to_study)

    #values should show starting with 1
    xvalues = range(len(values.values()[0]))
    xvalues.pop(0)
    xvalues.append(xvalues[-1] + 1)

    if not title:
        title = 'Nucleotide frequency per position'
    colors = {'A':'g', 'C':'b', 'T':'r', 'G':'k'}
    xlabel = 'Sequence positions'
    ylabel = 'Nucleotide frequency'

    draw_stacked_columns(values, colors, title=title, xlabel=xlabel,
                        ylabel=ylabel, fhand=fhand, xvalues=xvalues)
    return values
