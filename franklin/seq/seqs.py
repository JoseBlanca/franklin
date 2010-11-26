'''
Created on 2009 mar 27

@author: peio
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

from uuid import uuid4

import copy

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq as BioSeq
from Bio.Seq import UnknownSeq
from Bio.SeqFeature import SeqFeature as BioSeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Alphabet import (DNAAlphabet, ProteinAlphabet, Alphabet,
                          NucleotideAlphabet, SingleLetterAlphabet)

def copy_seq_with_quality(seqwithquality, seq=None, qual=None, name=None,
                          id_=None):
    '''Given a seqrecord it returns a new seqrecord with seq or qual changed.

    This is necessary because our SeqWithQuality is inmutable
    '''
    if seq is None:
        seq = seqwithquality.seq
    if id_ is  None:
        id_ = seqwithquality.id
    if name is None:
        name = seqwithquality.name

    #the letter annotations
    let_annot = {}
    for annot, value in seqwithquality.letter_annotations.items():
        if qual is not None and "phred_quality" == annot:
            let_annot[annot] = qual
        else:
            let_annot[annot] = value

    #the rest of parameters
    description = seqwithquality.description
    dbxrefs = seqwithquality.dbxrefs
    features = seqwithquality.features
    annotations = seqwithquality.annotations

    #the new sequence
    new_seq = SeqWithQuality(seq=seq, id=id_, name=name,
                             description=description, dbxrefs=dbxrefs,
                             features=features, annotations=annotations)

    #we restore the letter annotations (including the quality
    for annot, value in let_annot.items():
        new_seq.letter_annotations[annot] = value
    return new_seq

def get_seq_name(seq):
    'Given a sequence and its default name it returns its name'
    try:
        name = seq.name
    except AttributeError:
        name = None
    #the SeqRecord default
    if name == UNKNOWN_NAME:
        name = None

    if name is None:
        try:
            name = seq.id
        except AttributeError:
            name = None
    if name == UNKNOWN_ID:
        name = None

    if name is None:
        name = str(uuid4())
    return name

ALPHABETS = {'alphabet': Alphabet,
             'dnaalphabet': DNAAlphabet,
             'nucleotidealphabet':NucleotideAlphabet,
             'proteinalphabet':ProteinAlphabet,
             'singleletteralphabet':SingleLetterAlphabet}

_alphabet_name = lambda alpha: str(alpha).split('.')[-1].strip(')').strip('(')

ALPHABETS_REV = dict([(_alphabet_name(alpha), text) for text, alpha in ALPHABETS.items()])

UNKNOWN_NAME = "<unknown name>"
UNKNOWN_ID = "<unknown id>"
UNKNOWN_DESCRIPTION = "<unknown description>"

def _build_seq(struct):
    'Given a struct for a Seq it returns one Seq'
    if 'alphabet' in struct:
        alphabet = struct['alphabet']
        alphabet = ALPHABETS[alphabet]()
        return Seq(struct['seq'], alphabet=alphabet)
    else:
        return Seq(struct['seq'])

def create_seq_from_struct(struct):
    'It returns a SeqWithQuality'

    kwargs = {}

    kwargs['seq'] = _build_seq(struct['seq'])

    #the properties and default values
    properties = {'id': UNKNOWN_ID,
                  'name': UNKNOWN_NAME,
                  'description': UNKNOWN_DESCRIPTION,
                  'dbxrefs': None,
                  'annotations': None,
                  'letter_annotations':None}

    for prop, default in properties.items():
        value = struct[prop] if prop in struct else default
        kwargs[prop] = value

    #the features
    if 'features' not in struct:
        features = None
    else:
        features = []
        for feature in struct['features']:
            qualifiers = copy.deepcopy(feature['qualifiers'])
            type_ = feature['type']
            #The orf feature is special, it has two seqs in it
            if type_ == 'orf':
                if 'dna' in qualifiers:
                    qualifiers['dna'] = _build_seq(qualifiers['dna'])
                if 'pep' in qualifiers:
                    qualifiers['pep'] = _build_seq(qualifiers['pep'])
            feat = SeqFeature(location=FeatureLocation(int(feature['start']),
                                                       int(feature['end'])),
                              type=type_,
                              qualifiers=qualifiers)
            features.append(feat)
    kwargs['features'] = features

    return SeqWithQuality(**kwargs)

def fix_seq_struct_for_json(struct, alleles_to_string):
    'It returns a new alleles dict'
    for feature in struct['features']:
        #the snvs have tuples as keys and json doesn't like that
        if feature['type'] != 'snv':
            continue

        alleles = feature['qualifiers']['alleles']
        new_alleles = {}
        for allele_name, info in alleles.items():
            #we fix the name of the allele
            if alleles_to_string:
                allele_name = str(allele_name)
            else:
                allele_name = allele_name.strip().lstrip('(').rstrip(')')
                base, kind = allele_name.split(',')
                if base[0] == 'u':  #the unicode bit
                    base = base[1:]
                base = base.strip(base[0])
                kind = int(kind.strip())
                allele_name = base, kind
            new_alleles[allele_name] = info
        feature['qualifiers']['alleles'] = new_alleles

        #we fix the filter keys
        new_filters = {}
        if 'filters' in feature['qualifiers']:
            for filter_name, info in feature['qualifiers']['filters'].items():
                results = {}
                #info is a dict with a key for each param set
                for param_tuple, result in info.items():
                    if alleles_to_string:
                        param_tuple = repr(param_tuple)
                    else:
                        param_tuple = eval(param_tuple)
                    results[param_tuple] = result
                new_filters[filter_name] = results
            feature['qualifiers']['filters'] = new_filters
    return struct

def reverse_complement(seqrecord):
    'It return the reverse and complement sequence'
    try:
        qual = seqrecord.qual
    except AttributeError:
        qual = seqrecord.letter_annotations['phred_quality']

    seqwithqual =  SeqWithQuality(id=seqrecord.id, name=seqrecord.name,
                                  description=seqrecord.description,
                                  dbxrefs=seqrecord.dbxrefs,
                                  features=seqrecord.features,
                                  annotations=seqrecord.annotations,
                                  qual=qual, seq=seqrecord.seq)
    seqwithqual = seqwithqual.complement()

    return seqwithqual[::-1]

class SeqWithQuality(SeqRecord):
    '''A wrapper around Biopython's SeqRecord that adds a couple of convenience
    methods'''

    def __init__(self, seq, id=UNKNOWN_ID, name=UNKNOWN_NAME,
                 description=UNKNOWN_DESCRIPTION, dbxrefs=None,
                 features=None, annotations=None,
                 letter_annotations=None, qual=None):
        if id == UNKNOWN_ID and name != UNKNOWN_NAME:
            id = name
        #We don't want a Biopython Seq, we need our repr
        if not isinstance(seq, Seq) and not isinstance(seq, UnknownSeq):
            raise ValueError('seq should be a franklin Seq')
        SeqRecord.__init__(self, seq, id=id, name=name,
                           description=description, dbxrefs=dbxrefs,
                           features=features, annotations=annotations,
                           letter_annotations=letter_annotations)
        if qual is not None:
            self.qual = qual

    def _set_qual(self, qual):
        '''It stores the quality in the letter_annotations['phred_quality']'''
        self.letter_annotations["phred_quality"] = qual
    def _get_qual(self):
        '''It gets the quality from letter_annotations['phred_quality']'''
        if "phred_quality" not in self.letter_annotations:
            return None
        return self.letter_annotations["phred_quality"]
    qual = property(_get_qual, _set_qual)

    def complement(self):
        ''' it returns a new object with the complementary strand of the seq '''
        return self.__class__(seq=self.seq.complement(),
                              id=self.id + '_complemented',
                              name=self.name + '_complemented',
                              description=self.description,
                              dbxrefs=self.dbxrefs, features=self.features,
                              annotations=self.annotations,
                              letter_annotations=self.letter_annotations)

    def __add__(self, seq2):
        '''It returns a new object with both seq and qual joined '''
        #per letter annotations
        new_seq = self.__class__(name=self.name + '+' + seq2.name,
                                 id=self.id + '+' + seq2.id,
                                 seq=self.seq + seq2.seq)
        #the letter annotations, including quality
        for name, annot in self.letter_annotations.items():
            if name in seq2.letter_annotations:
                new_seq.letter_annotations[name] = annot + \
                                                   seq2.letter_annotations[name]
        return new_seq

    def __repr__(self):
        '''It writes the representation of the whole serecord,
        including annotations and feautures'''
        toprint = SeqRecord.__repr__(self)
        toprint = toprint[:-1]
        toprint += ', features=%s, ' % repr(self.features)
        toprint += 'annotations=%s, ' % repr(self.annotations)
        toprint += 'qual=%s,' % repr(self.qual)
        toprint += ")"
        return toprint

    def get_features(self, kind):
        'It yields the features that match the given kind'
        for feature in self.features:
            if feature.type == kind:
                yield feature

    def upper(self):
        'It returns the sequence upper cased'
        return self.__class__(seq=self.seq.upper(),
                              id=self.id,
                              name=self.name,
                              description=self.description,
                              dbxrefs=self.dbxrefs, features=self.features,
                              annotations=self.annotations,
                              letter_annotations=self.letter_annotations)

    def _from_seq_to_struct(self, seq):
        'Given a Seq it returns the struct'

        alphabet = ALPHABETS_REV[_alphabet_name(seq.alphabet)]

        struct = {'seq':str(seq)}
        if alphabet != 'alphabet': #the default one
            struct['alphabet'] = alphabet
        return struct

    def _get_struct(self):
        'It returns a structure with native python objects with all info'

        struct = {}

        struct['seq'] = self._from_seq_to_struct(self.seq)

        properties = {'id': UNKNOWN_ID,
                      'name': UNKNOWN_NAME,
                      'description': UNKNOWN_DESCRIPTION}
        for prop, default in properties.items():
            attr = getattr(self, prop)
            if attr != default:
                struct[prop] = attr

        #redundant information
        if 'id' in struct and struct['id'] == struct['name']:
            del struct['id']

        properties = ['dbxrefs', 'annotations', 'letter_annotations']
        for prop in properties:
            attr = getattr(self, prop)
            if attr:
                struct[prop] = attr

        features = []
        for feat in self.features:
            type_ = feat.type
            qualifiers = copy.deepcopy(feat.qualifiers)
            #The orf feature is special, it has two seqs in it
            if type_ == 'orf':
                if 'dna' in qualifiers:
                    qualifiers['dna'] = \
                                     self._from_seq_to_struct(qualifiers['dna'])
                if 'pep' in qualifiers:
                    qualifiers['pep'] = \
                                     self._from_seq_to_struct(qualifiers['pep'])
            loc = feat.location
            feat = {'start': loc.start.position,
                    'end':   loc.end.position,
                    'type':  type_,
                    'qualifiers': qualifiers}
            features.append(feat)
        struct['features'] = features

        return struct

    struct = property(_get_struct)

    def remove_annotations(self, kind):
        'It removes from seq annotations of any kind with any parameters'
        if kind == 'GOs':
            if  kind in self.annotations:
                del(self.annotations[kind])
        elif kind == 'orthologs':
            annotations = self.annotations
            for key in annotations.keys():
                if kind in key:
                    del(self.annotations[key])

        elif kind in ['microsatellite', 'orf', 'snv', 'intron']:
            new_features = []
            for feature in self.features:
                feat_kind = feature.type
                if feat_kind == kind:
                    continue
                new_features.append(feature)
            self.features = new_features
        elif kind == 'description':
            self.description = UNKNOWN_DESCRIPTION

class SeqFeature(BioSeqFeature):
    '''A wrapper around Biopython's SeqRecord that adds a couple of convenience
    methods'''
    def __init__(self, *args, **kwargs):
        BioSeqFeature.__init__(self, *args, **kwargs)

    def __repr__(self):
        'It prints representing the seqfeature'
        toprint = BioSeqFeature.__repr__(self)
        toprint = toprint[:-1]
        toprint += ', qualifiers=%s ' % repr(self.qualifiers)
        toprint += ")"
        return toprint


from Bio import Alphabet
from string import maketrans
from Bio.Data.IUPACData import (ambiguous_dna_complement,
                                ambiguous_rna_complement)
def _maketrans(complement_mapping) :
    """Makes a python string translation table (PRIVATE).

    Arguments:
     - complement_mapping - a dictionary such as ambiguous_dna_complement
       and ambiguous_rna_complement from Data.IUPACData.

    Returns a translation table (a string of length 256) for use with the
    python string's translate method to use in a (reverse) complement.

    Compatible with lower case and upper case sequences.

    For internal use only.
    """
    before = ''.join(complement_mapping.keys())
    after = ''.join(complement_mapping.values())
    before = before + before.lower()
    after = after + after.lower()
    return maketrans(before, after)

_dna_complement_table = _maketrans(ambiguous_dna_complement)
_rna_complement_table = _maketrans(ambiguous_rna_complement)

class Seq(BioSeq):
    'A biopython Seq with some extra functionality'
    def __eq__(self, seq):
        'It checks if the given seq is equal to this one'
        return str(self) == str(seq)

    def complement(self):
        """Returns the complement sequence. New Seq object.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_dna = Seq("CCCCCGATAG", IUPAC.unambiguous_dna)
        >>> my_dna
        Seq('CCCCCGATAG', IUPACUnambiguousDNA())
        >>> my_dna.complement()
        Seq('GGGGGCTATC', IUPACUnambiguousDNA())

        You can of course used mixed case sequences,

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_dna
        >>> my_dna = Seq("CCCCCgatA-GD", generic_dna)
        >>> my_dna
        Seq('CCCCCgatA-GD', DNAAlphabet())
        >>> my_dna.complement()
        Seq('GGGGGctaT-CH', DNAAlphabet())

        Note in the above example, ambiguous character D denotes
        G, A or T so its complement is H (for C, T or A).

        Trying to complement a protein sequence raises an exception.

        >>> my_protein = Seq("MAIVMGR", IUPAC.protein)
        >>> my_protein.complement()
        Traceback (most recent call last):
           ...
        ValueError: Proteins do not have complements!
        """
        base = Alphabet._get_base_alphabet(self.alphabet)
        if isinstance(base, Alphabet.ProteinAlphabet) :
            raise ValueError("Proteins do not have complements!")
        if isinstance(base, Alphabet.DNAAlphabet) :
            ttable = _dna_complement_table
        elif isinstance(base, Alphabet.RNAAlphabet) :
            ttable = _rna_complement_table
        elif ('U' in self._data or 'u' in self._data) \
        and ('T' in self._data or 't' in self._data):
            #TODO - Handle this cleanly?
            raise ValueError("Mixed RNA/DNA found")
        elif 'U' in self._data or 'u' in self._data:
            ttable = _rna_complement_table
        else:
            ttable = _dna_complement_table
        #Much faster on really long sequences than the previous loop based one.
        #thx to Michael Palmer, University of Waterloo
        return self.__class__(str(self).translate(ttable), self.alphabet)

    def upper(self):
        'It returns the uppercased sequence'
        return self.__class__(str(self).upper(), self.alphabet)

    def __getitem__(self, index) :                 # Seq API requirement
        #Note since Python 2.0, __getslice__ is deprecated
        #and __getitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        if isinstance(index, int) :
            #Return a single letter as a string
            return self._data[index]
        else :
            #Return the (sub)sequence as another Seq object
            return self.__class__(self._data[index], self.alphabet)

    def __repr__(self):
        """Returns the representation of the sequence."""
        return "%s(%s, %s)" % (self.__class__.__name__,
                               #repr(self.data),
                               repr(self.tostring()),
                               repr(self.alphabet))
