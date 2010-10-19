'''
Created on Oct 4, 2010

@author: peio
'''

class CodonPosition(object):
    'It represents a position in a codon'
    def __init__(self, codon, position):
        '''It inits the object.

        codon - codon number (int)
        position - position in the codon (0, 1 or 2)
        '''
        self.num = self._to_num((codon, position))

    def _get_codon(self):
        'It returns the codon number'
        codon = abs(self.num) // 3
        if self.num < 0:
            codon = -codon
        return codon

    def _get_position(self):
        'It returns the position in the codon'
        position = abs(self.num) % 3
        if self.num < 0:
            position = - position
        return position
    codon = property(_get_codon)
    position = property(_get_position)

    @staticmethod
    def _to_num(obj):
        'It converts to number'
        if isinstance(obj, int):
            return obj
        elif isinstance(obj, tuple):
            return obj[0] * 3 + obj[1]
        else:
            return obj.num

    def __eq__(self, codon):
        'It returns True if both codons are equal'
        return self.num == self._to_num(codon)

    def __ne__(self, codon):
        'It returns True if both codons are not equal'
        return self.num != self._to_num(codon)

    def __lt__(self, codon):
        'It returns True if both codons are equal'
        return self.num < self._to_num(codon)

    def __le__(self, codon):
        'It returns True if both codons are equal'
        return self.num <= self._to_num(codon)

    def __gt__(self, codon):
        'It returns True if both codons are equal'
        return self.num > self._to_num(codon)

    def __ge__(self, codon):
        'It returns True if both codons are equal'
        return self.num >= self._to_num(codon)

    def __add__(self, obj):
        'It adds an int or a CodonPosition'
        return self.__class__(0, self._to_num(self) + self._to_num(obj))

    def __radd__(self, obj):
        'It adds an int or a CodonPosition'
        return self.__class__(0, self._to_num(self) + self._to_num(obj))

    def __sub__(self, obj):
        'The substraction'
        return self.__class__(0, self._to_num(self) - self._to_num(obj))

    def __rsub__(self, obj):
        'The substraction'
        return self.__class__(0, self._to_num(self) - self._to_num(obj))

    def __neg__(self):
        'Unary negation'
        return self.__class__(0, self.num)

    def __abs__(self):
        'absolute value'
        return self.__class__(0, abs(self.num))

    def __nonzero__(self):
        'False value for object'
        return bool(self.num)

    def __str__(self):
        'The string representation'
        return '%i %i' % (self.codon, self.position)

    def __repr__(self):
        'The string representation'
        return 'CodonPosition(%i %i)' % (self.codon, self.position)

class CoordSystem(object):
    'It represents the different coordinate systems in a gene'
    def __init__(self, relations):
        '''The init

        It requires a dict with the fragments. Each key is the name of the gene,
        cDNA or protein. The values are lists of tuples. Each tuple has the
        beginning and the end of the fragment. e.g:
          {'genomic': [(2, 5),  (8, 11)],
           'cdna':    [(0, 3),  (4, 7) ],
           'prot':    [((0, 0), (0, 1)),
                       ((0, 2), (1, 2))]}

        #genomic           111
        #        0123456789012
        #
        #cdna      0123  4567
        #
        #prot        00  0111
        #            01  2012

        The protein positions are represented by a tuple, the codon and the
        position in the codon.
        A fragment can also be None if it is not defined for this molecule.
        '''
        self._fragments = None
        self._process_relations(relations)

    def _add_relation(self, fragments, relation):
        'It adds one relation to the fragments'

        #definitions
        # segment. a tuple of limits that define a region
        # fragment. a segment with the information for all molecules
        # limit. an integer or CodonPosition that defines a segment

        molecules = relation.keys()
        if len(molecules) != 2:
            raise 'relation with more than two molecules: ' + str(relation)

        mol0, mol1 = molecules

        #are the segments defined by the two molecules in the relations equal?
        num_segments = len(relation.values()[0])
        for index in range(num_segments):
            seg0 = relation[mol0][index]
            seg1 = relation[mol1][index]
            if abs(seg0[1] - seg0[0]) != abs(seg1[1] - seg1[0]):
                msg = 'Relation with inequal lengths: ' + str(relation)
                raise ValueError(msg)

        if not fragments:
            fragments = self._create_fragments_from_relation(relation, mol0,
                                                             mol1)
            return fragments

        #which molecule is common between the new relation and the old
        #fragments?
        molecules_in_frags = set()
        for frag in fragments:
            molecules_in_frags.update(frag.keys())
        if mol0 in molecules_in_frags and mol1 not in molecules_in_frags:
            ref_mol = mol0
            new_mol = mol1
        elif mol1 in molecules_in_frags and mol0 not in molecules_in_frags:
            ref_mol = mol1
            new_mol = mol0
        elif mol1 in molecules_in_frags and mol0 in molecules_in_frags:
            msg = 'Both molecules were in the previous relations: %s' % \
                                                                   str(relation)
            raise ValueError(msg)
        else:
            msg = 'No molecule were in the previous relations: %s' % \
                                                                   str(relation)
            raise ValueError(msg)

        #which are the old segments
        old_segments = [frag[ref_mol] for frag in fragments if ref_mol in frag]
        rel_segments = relation[ref_mol]
        new_segments = self._intersec_segments(old_segments, rel_segments)

        #we have to split the old fragments following the new segments
        new_fragments = []
        for frag in fragments:
            if ref_mol in frag:
                frag_segment = frag[ref_mol]
                #which of the new segments are inside this fragment
                new_frag_segments = self._get_segments_in_segment(new_segments,
                                                                  frag_segment)
                #we divide the fragment into the new segments
                divided_frags = self._split_fragment(frag, new_frag_segments,
                                                     ref_mol)
                new_fragments.extend(divided_frags)
            else:
                new_fragments.append(frag)

        #now we can add the new relation to these new fragments
        self._add_relation_to_fragments(relation, new_mol, ref_mol,
                                        new_fragments)

        return new_fragments

    @staticmethod
    def _get_segments_in_segment(segments, segment):
        'It returns the segments that are inside the given segment'

        segment0, segment1 = min(segment), max(segment)

        inside_segments = []
        for seg in segments:
            seg0 = min(seg)
            if segment0 <= seg0 <= segment1:
                inside_segments.append(seg)
        return inside_segments

    @staticmethod
    def _merge_overlaping_segments(segments):
        'It merges the segments'
        segs = segments[:]
        segs.sort(key=lambda x: x[0])

        merged_segs = [segs[0]]
        for seg in segs:
            seg0 = merged_segs[-1]
            seg1 = seg
            #do this one and the next overlap
            if seg0[1] >= seg1[0]:
                new_seg = min(seg0[0], seg1[0]), max(seg0[1], seg1[1])
                merged_segs[-1] = new_seg
            else:
                merged_segs.append(seg1)
        return merged_segs

    def _intersec_segments(self, segments1, segments2):
        '''Given two lists of segments it returns the intersected list

        A segment is a pair of numbers that define a start and an end.
            <----->
        This function does:
             fragments1 <-----------------><--------------->
             fragments2     <-----><----------- ><--------->
             result     <--><-----><------><----><--------->
        '''

        reverse_segments = lambda x: [(e[1], e[0]) for e in reversed(x)]
        #the following code only works with forward segments
        #so we reverse the reverse fragments
        forward = True
        if segments1[0][1] - segments1[0][0] < 0:
            segments1 = reverse_segments(segments1)
            forward = False
        if segments2[0][1] - segments2[0][0] < 0:
            segments2 = reverse_segments(segments2)

        #which are the regions covered by these segments
        segs = segments1[:]
        segs.extend(segments2)
        covered_regions = self._merge_overlaping_segments(segs)

        #we collect the fragment limits
        limits = set()
        for segments in (segments1, segments2):
            for segment in segments:
                limits.add((segment[0], '<'))
                limits.add((segment[1], '>'))
        # Now we sort the limits
        limits = list(limits)
        limits.sort(key=lambda x: x[1])
        limits.sort(key=lambda x: x[0])

        #the open and close limits should alternate
        #we add the required limits to comply
        limit_index = 0
        while limit_index < len(limits):
            if (limits[limit_index][1] == '<' and
                limits[limit_index + 1][1] == '>'):
                #all ok, we move to the next one
                limit_index += 2
            elif (limits[limit_index][1] == '<' and
                limits[limit_index + 1][1] == '<'):
                #we have to add a new closing limit before the next opening one
                limits.insert(limit_index + 1,
                              (limits[limit_index + 1][0] - 1, '>'))
                limit_index += 2
            elif limits[limit_index][1] == '>':
                limits.insert(limit_index ,
                              (limits[limit_index - 1][0] + 1, '<'))

        #is there any gaps
        # <  ><  ><  >  GAP    <   ><   >
        limit_index = 0
        while limit_index < len(limits):
            pair0 = limits[limit_index], limits[limit_index + 1]
            try:
                pair1 = limits[limit_index + 2], limits[limit_index + 3]
            except IndexError:
                break
            #is there a gap here
            if pair1[0][0] - pair0[1][0] <= 1:
                limit_index += 2
                continue
            #we create a new segment to cover the gap
            new_segment = pair0[1][0] + 1, pair1[0][0] - 1
            #does the new segment overlaps with any of the covered regions?
            overlaps = False
            for reg in covered_regions:
                if reg[0] <= new_segment[0] <= reg[1]:
                    overlaps = True
            if overlaps:
                #we add the new segment to cover the gap
                limits.insert(limit_index + 2, (new_segment[0], '<'))
                limits.insert(limit_index + 3, (new_segment[1], '>'))
                limit_index += 2
            limit_index += 2

        #now we can create the new segments
        segments = []
        limit_index = 0
        while limit_index < len(limits):
            open_limit = limits[limit_index]
            close_limit = limits[limit_index + 1]
            if open_limit[1] != '<':
                raise RuntimeError('Opening limit is closing limit')
            if close_limit[1] != '>':
                raise RuntimeError('Closing limit is opening limit')
            segments.append((open_limit[0], close_limit[0]))
            limit_index += 2

        if not forward:
            segments = reverse_segments(segments)

        return segments

    @staticmethod
    def _create_fragments_from_relation(relation, molecule_0, molecule_1):
        'It creates the fragments given one relation'
        fragments = []
        for mol_frag0, mol_frag1 in zip(relation[molecule_0],
                                        relation[molecule_1]):
            fragment = {}
            fragment[molecule_0] = mol_frag0
            fragment[molecule_1] = mol_frag1
            fragments.append(fragment)
        return fragments

    def _get_frag_limits_and_segments(self, fragments, molecule):
        '''Given some fragments it returns the list of limits and the fragment
        limits

        The limits are just ints and the fragment limits are tuples
        '''
        limits = set()
        frag_limits = []
        for frag in fragments:
            if molecule not in frag:
                continue
            limits.update(frag[molecule])
            frag_limits.append(frag[molecule])
        frag_limits.sort(key=lambda x: x[0])
        limits = sorted(list(limits))
        return limits, frag_limits

    def _add_relation_to_fragments(self, relation, new_molecule, ref_molecule,
                                   fragments):
        '''Given a relation and the fragments it adds the relation to those
        fragments.

        The relation relates two molecules, only the one specified by the
        given new_molecule will be added.
        '''

        #which would we the fragments if these were the only relation?
        rel_frags = self._create_fragments_from_relation(relation, new_molecule,
                                                         ref_molecule)
        #we have to split the relation fragments into the fragment limits
        #defined by the fragments
        frag_limits = self._get_frag_limits_and_segments(fragments,
                                                         ref_molecule)[1]
        new_rel_frags = []
        for rel_frag in rel_frags:
            new_rel_frags.extend(self._split_fragment(rel_frag, frag_limits,
                                         ref_molecule))
        rel_frags = new_rel_frags

        #now we add these fragments to the fragments, but only the new_molecule
        for rel_frag in rel_frags:
            #which of the old fragments correspond to this frag derived from the
            #relation
            old_frag = self._get_fragment_with_limit(fragments, ref_molecule,
                                                     rel_frag[ref_molecule])
            old_frag[new_molecule] = rel_frag[new_molecule]

    def _split_fragment(self, fragment, new_segments, new_segment_molecule):
        'Given a fragment a some new limits it splits the fragment'

        #the limits for the reference molecule
        ref_segment = fragment[new_segment_molecule]
        new_fragments = []
        for new_segment in new_segments:
            new_fragment = {}
            for molecule, limit in fragment.items():
                #if the new_segment is outside the boundaries of the reference
                #limits for this fragment we do not add it here.
                new_seg0, new_seg1 = min(new_segment), max(new_segment)
                ref_seg0, ref_seg1 = min(ref_segment), max(ref_segment)
                if ((new_seg1 < ref_seg0) or
                    (new_seg0 > ref_seg1)):
                    continue
                if molecule == new_segment_molecule:
                    pos0, pos1 = new_segment
                else:
                    pos0 = self._transform_pos(ref_segment, limit,
                                               new_segment[0])
                    pos1 = self._transform_pos(ref_segment, limit,
                                               new_segment[1])
                new_fragment[molecule] = (pos0, pos1)
            if new_fragment:
                new_fragments.append(new_fragment)
        return new_fragments

    @staticmethod
    def _transform_pos(limits1, limits2, position1):
        '''Given a position in the coord of the molecule1 it returns the
        position in the coord of the molecule2'''

        start1 = limits1[0]
        end1   = limits1[1]
        start2 = limits2[0]
        end2   = limits2[1]

        #is forward or reverse?
        if end1 - start1 == end2 - start2:
            forward = True
        else:
            forward = False
        if forward:
            pos2 = position1 - start1 + start2
        else:
            pos2 = start2 + start1 - position1

        #do we return an int or a CodonPosition
        if isinstance(limits2[0], int) and not isinstance(pos2, int):
            pos2 = pos2.num
        elif not isinstance(limits2[0], int) and isinstance(pos2, int):
            pos2 = CodonPosition(codon=0, position=pos2)
        return pos2

    @staticmethod
    def _get_fragment_with_limit(fragments, molecule, limits):
        'It returns the fragments with the given limits for the given molecule'
        limits = sorted(limits)
        for frag in fragments:
            if molecule not in frag:
                continue
            frag_limits = frag[molecule]
            frag_limits = sorted(frag_limits)
            if limits == frag_limits:
                return frag
        msg = 'fragment not found for molecule %s and limits %s' % (molecule,
                                                                    str(limits))
        raise RuntimeError(msg)

    def _process_relations(self, relations):
        'It transforms the fragments into the structure used by the class'

        frags = []
        for relation in relations:
            frags = self._add_relation(frags, relation)
        self._fragments = frags

    def get_fragment(self, molecule, position):
        'It returns the fragment that covers the given position'

        for frag in self._fragments:
            if molecule not in frag:
                continue
            range_ = frag[molecule]
            if range_[0] <= position <= range_[1]:
                return frag
        raise RuntimeError('fragment not found for: %s, %s' % (molecule,
                                                               str(position)))

    def get_molecule_min_position(self, molecule):
        'It returns the start of the coord system of the molecule'

        mol_min = None
        for frag in self._fragments:
            if molecule not in frag:
                continue
            this_min = min(frag[molecule])
            if mol_min is None or mol_min > this_min:
                mol_min = this_min
        return mol_min

    def transform(self, from_mol, to_mol, position):
        'It transforms the coord from one molecule to the other'

        #which is the fragment for this molecule and position
        fragment = self.get_fragment(molecule=from_mol, position=position)

        return self._transform_pos(fragment[from_mol], fragment[to_mol],
                                   position)
