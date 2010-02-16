'''An fpc physical map representation with gff export capabilities.

This module has been coded looking at the bioperl equivalent module
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

import re
from franklin.gff import gff_parser

def fpcgff2_parser(fhand):
    'It parses a gff2_annotations file and yields converted gff3 features'

    for feature in gff_parser(fhand, 2):
        annotations = feature['attributes']
        name        = annotations['Name']
        type_ = feature['type']

        feature['id']      = name
        feature['name']    = name

        if type_ in ('Chromosome','contig', 'marker'):
            pass

        elif type_ == 'BAC':
            contig_name = 'ctg%d' % int(annotations['Contig_hit'])

            if 'Marker_hit' in annotations:
                feature['marker_hit'] = annotations['Marker_hit']

        else:
            raise ValueError('Unknown feature type: %s' % type_)

        del feature['attributes']

        yield feature

class FPCMap(object):
    '''It parses and hold the information from an FPC map'''

    def __init__(self, fpc_file):
        '''It parses and fpc file'''
        self._fhand = fpc_file
        self.name = None
        self.version = None
        self.clones = {}
        self.contigs = {}
        self.markers = {}
        self._parse_fpc()

    def _parse_fpc(self):
        'It parses the fpc file'
        #the header
        version_re = re.compile('^((\d+\.\d+)).*')
        clone_lpos_re = re.compile('Map "ctg(\d+)" Ends Left ([-\d]+)')
        clone_rpos_re = re.compile('Map "ctg\d+" Ends Right ([-\d]+)')
        clone_match_re = re.compile('([a-zA-Z]+)_match_to_\w+\s+"(.+)"')
        clone_positive_re = re.compile('Positive_(\w+)\s+"(.+)"')
        for line in self._fhand:
            line = line.strip()
            if not line:
                break
            line = line.lstrip('/ ')
            #project name
            if 'project' in line and 'fpc' in line:
                self.name = line.split()[-1]
            version_match = version_re.match(line)
            if version_match:
                self.version = line.split()[0]

        #clone data
        clones = self.clones
        contigs = self.contigs
        markers = self.markers
        for line in self._fhand:
            if line.startswith('Markerdata'):
                break
            line = line.strip()
            if not line:
                continue
            type_, name = line.split(':')
            clone = {}
            clone['type'] = type_.strip()
            clone['name'] = name.strip().strip('"')
            #clone attributes
            fpc_remark = ''
            remark = ''
            for line in self._fhand:
                line = line.strip()
                if not line:
                    break
                lpos_match = clone_lpos_re.match(line)
                if lpos_match:
                    contig_number = lpos_match.group(1)
                    contig_start = lpos_match.group(2)
                    if 'range' not in clone:
                        clone['range'] = {}
                    clone['contig_number'] = contig_number
                    clone['range']['start'] = contig_start
                    continue
                rpos_match = clone_rpos_re.match(line)
                if rpos_match:
                    contig_end = rpos_match.group(1)
                    if 'range' not in clone:
                        clone['range'] = {}
                    clone['range']['end'] = contig_end
                    continue
                match_match = clone_match_re.match(line)
                if match_match:
                    match_type = 'match_' + match_match.group(1).lower()
                    matched_object = match_match.group(2)
                    if match_type not in clone:
                        clone[match_type] = []
                    clone[match_type].append(matched_object)
                    continue
                positive_match = clone_positive_re.match(line)
                if positive_match:
                    if 'markers' not in clone:
                        clone['markers'] = []
                    marker_type = positive_match.group(1)
                    marker = positive_match.group(2)
                    clone['markers'].append(marker)
                    if marker not in markers:
                        markers[marker] = {}
                        markers[marker]['contigs'] = {}
                        markers[marker]['clones'] = []
                    markers[marker]['type'] = marker_type
                    markers[marker]['contigs'][clone['contig_number']] = True
                    markers[marker]['clones'].append(clone['name'])
                elif line.startswith('Gel_number'):
                    gel_number = line.split()[1]
                    clone['gel'] = gel_number
                elif line.startswith('Remark'):
                    remark += line.split(' ', 1)[1].strip('"')
                    remark += '\n'
                    if 'Chr' in remark:
                        raise NotImplemented('Fixme')
                elif line.startswith('Fp_number'):
                    fp_number = line.split()[1]
                    clone['fp_number'] = fp_number
                elif line.startswith('Shotgun'):
                    seq_type, seq_status = line.split()[1:]
                    clone['sequence_type'] = seq_type
                    clone['sequence_status'] = seq_status
                elif line.startswith('fpc_remark'):
                    fpc_remark += line.split(' ', 1)[1].strip('"')
                    fpc_remark += '\n'
            clone['remark'] = remark
            clone['fpc_remark'] = fpc_remark
            clones[clone['name']] = clone
            if 'contig_number' in clone:
                contig_number = clone['contig_number']
                if contig_number not in contigs:
                    contigs[contig_number] = {}
                    contigs[contig_number]['clones'] = []
                    contigs[contig_number]['range'] = {}
                    contigs[contig_number]['range']['start'] = None
                    contigs[contig_number]['range']['end'] = None
                #this clone is in the contig
                contigs[contig_number]['clones'].append(clone['name'])
                #is the contig larger due to this clone?
                start = clone['range']['start']
                end = clone['range']['end']
                if (contigs[contig_number]['range']['start'] is None or
                    contigs[contig_number]['range']['start'] > start):
                    contigs[contig_number]['range']['start'] = start
                if (contigs[contig_number]['range']['end'] is None or
                    contigs[contig_number]['range']['end'] < end):
                    contigs[contig_number]['range']['end'] = end
            #the markers in the contig
            if 'markers' in clone:
                if 'markers' not in contigs[contig_number]:
                    contigs[contig_number]['markers'] = {}
                for marker in clone['markers']:
                    contigs[contig_number]['markers'][marker] = True

        #marker data
        re_group = re.compile('(\d+|\w)(.*)')
        #re_anchor = re.compile('Anchor_pos\s+([\d.]+)\s+(F|P)?')
        for line in self._fhand:
            if line.startswith('Contigdata'):
                break
            line = line.strip()
            if not line:
                continue
            type_, name = line.split(':')
            name = name.strip().strip('"')
            type_ = type_.split('_', 1)[-1]
            marker = {}
            marker['type'] = type_
            marker['group'] = 0
            marker['global'] = 0
            marker['anchor'] = 0
            remark = ''
            for line in self._fhand:
                line = line.strip()
                if not line:
                    break
                if line.startswith('Global_position'):
                    raise NotImplemented('Fixme')
                elif line.startswith('Anchor_bin'):
                    group = line.split()[1].strip('"')
                    match = re_group.match(group)
                    marker['group'] = match.group(1)
                    marker['subgroup'] = match.group(2)
                elif line.startswith('Anchor_pos'):
                    match = re_group.match(line)
                    marker['global'] = match.group(1)
                    marker['anchor'] = 1
                    if match.group(2) == 'F':
                        marker['framework'] = 1
                    else:
                        marker['framework'] = 0
                elif line.startswith('anchor'):
                    marker['anchor'] = 1
                elif line.startswith('Remark'):
                    remark += line.split(' ', 1)[1].strip('"')
                    remark += '\n'
            markers[name] = marker

        #contig data
        re_contig_def = re.compile('^Ctg(\d+)')
        re_contig_remark = re.compile('#\w*(.*)\w*$/')
        re_chr_remark = re.compile('Chr_remark\s+"(-|\+|Chr(\d+))\s+(.+)"$')
        for line in self._fhand:
            line = line.strip()
            if not line:
                continue
            if line.startswith('Ctg'):
                match = re_contig_def.match(line)
                contig = None
                name = match.group(1)
                if name in contigs:
                    contig = contigs[name]
                else:
                    contig = {}
                    contigs[name] = contig
                contig['group'] = 0
                contig['anchor'] = 0
                contig['position'] = 0
                match = re_contig_remark.search(line)
                if match:
                    contig['group'] = match.group(1)
                    contig['anchor'] = 1
            elif line.startswith('Chr_remark'):
                match = re_chr_remark.match(line)
                if match is None:
                    continue
                contig['anchor'] = 1
                if match.group(3):
                    contig['chr_remark'] = match.group(3)
                if match.group(2):
                    contig['group'] = match.group(2)
                else:
                    contig['group'] = '?'
            elif line.startswith('User_remark'):
                remark = line.split()[1].strip('"')
                contig['usr_remark'] = remark
            elif line.startswith('Trace_remark'):
                remark = line.split()[1].strip('"')
                contig['trace_remark'] = remark

    #relevant SO terms
    #BAC
    #A clone is composed by a vector replicon and a clone_insert. The SO BAC is
    #the vector_replicon. The clone_insert (SO:0000753), in this case, is a
    #cloned_genomic_insert (SO:0000914) and more specifically a
    #BAC_cloned_genomic_insert (SO:0000992).

    #Physical map
    #The map is a fragment_assembly (SO:0001249), although there is a more
    #specific fingerprint_map (SO:0001250)

    #[Term]
    #id: SO:0000151
    #name: clone
    #def: "A piece of DNA that has been inserted in a vector so that it can be
    #      propagated in a host bacterium or some other organism." [SO:ke]
    #subset: SOFA
    #xref: http:http\://en.wikipedia.org/wiki/Clone_(genetics) "wiki"
    #is_a: SO:0000695 ! reagent

    #[Term]
    #id: SO:0000153
    #name: BAC
    #def: "Bacterial Artificial Chromosome, a cloning vector that can be
    #      propagated as mini-chromosomes in a bacterial host." [SO:ma]
    #comment: This term is mapped to MGED. Do not obsolete without consulting
    #         MGED ontology.
    #synonym: "bacterial artificial chromosome" EXACT []
    #is_a: SO:0000440 ! vector_replicon
    #is a SO:0000440(vector_replicon)

    #[Term]
    #id: SO:0000753
    #name: clone_insert
    #def: "The region of sequence that has been inserted and is being
    #      propogated by the clone." [SO:ke]
    #synonym: "clone insert" EXACT []
    #is_a: SO:0000695 ! reagent
    #relationship: part_of SO:0000151 ! clone

    #[Term]
    #id: SO:0000992
    #name: BAC_cloned_genomic_insert
    #comment: Requested by Andy Schroder - Flybase Harvard, Nov 2006.
    #synonym: "BAC cloned genomic insert" EXACT []
    #is_a: SO:0000914 ! implied link automatically realized
    #intersection_of: SO:0000914 ! cloned_genomic_insert
    #intersection_of: derives_from SO:0000153 ! BAC
    #relationship: derives_from SO:0000153 ! implied link automatically realized

    #[Term]
    #id: SO:0001249
    #name: fragment_assembly
    #def: "A fragment assembly is a genome assembly that orders overlapping
    #      fragments of the genome based on landmark sequences. The base pair
    #      distance between the landmarks is known allowing additivity of
    #      lengths." [SO:ke]
    #synonym: "fragment assembly" EXACT []
    #synonym: "physical map" EXACT []
    #is_a: SO:0001248 ! assembly

    #[Term]
    #id: SO:0001250
    #name: fingerprint_map
    #def: "A fingerprint_map is a physical map composed of restriction
    #     fragments." [SO:ke]
    #synonym: "BACmap" EXACT []
    #synonym: "fingerprint map" EXACT []
    #synonym: "FPC" EXACT []
    #synonym: "FPCmap" EXACT []
    #synonym: "restriction map" EXACT []
    #is_a: SO:0001249 ! fragment_assembly
    #relationship: has_part SO:0000412 ! restriction_fragment

    #[Term]
    #id: SO:0001251
    #name: STS_map
    #def: "An STS map is a physical map organized by the unique STS landmarks."
    #synonym: "STS map" EXACT []
    #is_a: SO:0001249 ! fragment_assembly
    #relationship: has_part SO:0000331 ! STS
