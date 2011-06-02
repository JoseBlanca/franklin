'''
Created on 06/05/2011

@author: jose
'''
import random

from franklin.seq.seqs import SeqWithQuality, Seq

def create_random_seq(length, gc=50):
    'It returns a random sequence'
    nucl_choice = {'at': {0:'A', 1:'T'}, 'gc': {0:'G', 1:'C'}}
    gc = gc / 100.0    
    nucls = []
    for index in range(length):
        gc_choice = random.uniform(0, 1)
        gc_choice = 'gc' if gc_choice <= gc else 'at'
        at_choice = round(random.uniform(0, 1))
        nucl = nucl_choice[gc_choice][at_choice]        
        nucls.append(nucl)        
    assert len(nucls) == length
    return ''.join(nucls)

def create_random_seqwithquality(length, gc=50, qual_range=40):
    'It returns a random seqwithquality'
    if isinstance(qual_range, int):
        qual_range = [qual_range, qual_range]
        quals = []
    for index in range(length):
        qual = random.randint(qual_range[0], qual_range[1])
        quals.append(qual)
    seq = create_random_seq(length, gc)
    name = list('holacaracola')
    random.shuffle(name)
    name = ''.join(name)
    return SeqWithQuality(Seq(seq), qual=quals, name=name)

