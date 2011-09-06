'''
Created on 18/07/2011

@author: dani
'''

from franklin.seq.readers import seqs_in_file
from franklin.snv.snv_annotation import (calculate_maf_frequency,
                                         calculate_heterozygosity,
                                         calculate_pic)
from os.path import join
from optparse import OptionParser
from operator import itemgetter

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-s', '--inseqfile', dest='inseqfile',
                      help='input sequence file, (pickle)')

    parser.add_option('-d', '--dir_out', dest='dir_out',
                      help='Output directory')

    parser.add_option('-g', '--group_kind', dest='group_kind',
        help='Group kind to use (library(LB), platform(PL), sample(SM))')

    parser.add_option('-W', '--window_width', dest='window_width',
        help='Width of the window')

    parser.add_option('-w', '--window_step', dest='window_step',
        help='Displacement of the window')

    parser.add_option('-v', '--value_kind', dest='value_kind',
        help='Kind of value to be drawn: het, pic or maf ')

    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.inseqfile:
        in_fpath = options.inseqfile
    else:
        parser.error('Input seq file mandatory')

    if options.value_kind:
        value_kind = options.value_kind
        if value_kind not in ('het', 'pic', 'maf'):
            parser.error('Value kind must be het, pic or maf')
    else:
        parser.error('Value kind mandatory')

    if options.dir_out:
        dir_out = options.dir_out
    else:
        dir_out = '.'

    if options.group_kind:
        group_kind = options.group_kind
        if group_kind not in ('LB', 'PL', 'SM'):
            parser.error('Group kind must be LB, PL or SM')
    else:
        group_kind = None

    if options.window_width:
        window_width = int(options.window_width)
    else:
        window_width = None

    if options.window_step:
        window_step = int(options.window_step)
    else:
        window_step = None

    return (in_fpath, dir_out, group_kind, window_width, window_step,
            value_kind)

def get_groups(fpath):
    'It gets the values of every kind of group'
    groups = {'LB':[], 'PL':[], 'SM':[]}

    for seq in seqs_in_file(open(fpath)):
        for snv in seq.get_features(kind='snv'):
            read_groups = snv.qualifiers['read_groups']
            for tags in read_groups.values():
                for group_kind, group in tags.items():
                    if group not in groups[group_kind]:
                        groups[group_kind].append(group)
    return groups

def calculate_hets_group(seqs, groups=None, group_kind=None, ploidy=2):
    'It calculates the snv heterozygosity of a given group'
    het_profile = {}
    for seq in seqs:
        for snv in seq.get_features('snv'):
            het = calculate_heterozygosity(snv, ploidy,
                                           group_kind=group_kind,
                                           groups=groups)
            if het is not None:
                location = snv.location.start.position
                seq_name = seq.name
                if seq_name not in het_profile:
                    het_profile[seq_name] = []
                het_profile[seq_name].append((location, het))
    return het_profile

def calculate_pics_group(seqs, groups=None, group_kind=None):
    'It calculates the snv heterozygosity of a given group'
    pic_profile = {}
    for seq in seqs:
        for snv in seq.get_features('snv'):
            pic = calculate_pic(snv, group_kind=group_kind, groups=groups)
            if pic is not None:
                location = snv.location.start.position
                seq_name = seq.name
                if seq_name not in pic_profile:
                    pic_profile[seq_name] = []
                pic_profile[seq_name].append((location, pic))
    return pic_profile

def calculate_mafs_group(seqs, groups=None, group_kind=None):
    'It calculates the snv heterozygosity of a given group'
    maf_profile = {}
    for seq in seqs:
        for snv in seq.get_features('snv'):
            maf = calculate_maf_frequency(snv, group_kind=group_kind,
                                          groups=groups)
            if maf is not None:
                location = snv.location.start.position
                seq_name = seq.name
                if seq_name not in maf_profile:
                    maf_profile[seq_name] = []
                maf_profile[seq_name].append((location, maf))
    return maf_profile

def apply_window(profile, window_width, window_step):
    'It modifies the profile, calculating values within a given moving window'
    if window_width % 2 == 0:
        raise ValueError('Window width must be an odd number, not even number')

    window_arm = (window_width - 1 ) / 2

    new_profile = {}
    for sequence in profile:
        sequence_profile = profile[sequence]
        if len(sequence_profile) > 1:
            sorted(sequence_profile, key=itemgetter(0))

        if sequence_profile[0][0] < (window_width - 1):
            window_center = window_arm
        else:
            window_center = sequence_profile[0][0] - window_arm

        left_margin = window_center - window_arm
        right_margin = window_center + window_arm

        n = 0
        snvs_in_window = []
        snvs_in_next_window = []
        new_sequence_profile = []
        while (right_margin <= sequence_profile[-1][0] and
               n < len(sequence_profile)):
            if (sequence_profile[n][0] >= left_margin and
                sequence_profile[n][0] <= right_margin):
                snvs_in_window.append(sequence_profile[n][1])
                #is the snv also in the next window?:
                if sequence_profile[n][0] >= (left_margin + window_step):
                    snvs_in_next_window.append(sequence_profile[n][1])
                n += 1
            else:
                if len(snvs_in_window) >= 3:
                    average = sum(snvs_in_window)/len(snvs_in_window)
                    new_sequence_profile.append((window_center, average))

                right_margin += window_step
                window_center = right_margin - window_arm
                left_margin = window_center - window_arm

                snvs_in_window = []
                for snv in snvs_in_next_window:
                    snvs_in_window.append(snv)
                snvs_in_next_window = []

        new_profile[sequence] = new_sequence_profile
    return new_profile

def write_wig(dir_out, profile, group, window_width, window_step):
    'It writes the wig file given a profile'
    wig_fhand = open(join(dir_out, group+'.wig'), 'w')

    for sequence in profile:
        header = ('variableStep'+'\t'+'chrom='+sequence+'\t'+'span='+
                  str(window_width)+'\n')
        wig_fhand.write(header)
        sequence_profile = profile[sequence]
        for snv in sequence_profile:
            if window_width:
                window_arm = (window_width - 1 ) / 2
                # 1-indexed
                left_margin = snv[0] - window_arm + 1
                wig_fhand.write(str(left_margin)+'\t'+str(snv[1])+'\n')
            else:
                wig_fhand.write(str(snv[0])+'\t'+str(snv[1])+'\n')
    wig_fhand.flush()

def draw_sequence_distribution(seq_path, dir_out, group_kind,
                               window_width, window_step,
                               value_kind):
    '''It creates a wig file, given a value kind, with the distribution of that
    value along a given sequence'''
    groups = get_groups(seq_path)
    for group in groups[group_kind]:
        seqs = seqs_in_file(open(seq_path, 'r'))

        if value_kind == 'het':
            profile = calculate_hets_group(seqs, groups=[group],
                                           group_kind=group_kind)
        if value_kind == 'pic':
            profile = calculate_pics_group(seqs, groups=[group],
                                           group_kind=group_kind)
        if value_kind == 'maf':
            profile = calculate_mafs_group(seqs, groups=[group],
                                           group_kind=group_kind)

        if profile and window_width and window_step:
            new_profile = apply_window(profile,
                                       window_width=window_width,
                                       window_step=window_step)

            write_wig(dir_out, new_profile, group, window_width,
                      window_step)
        else:
            write_wig(dir_out, profile, group)

def main():
    'The main part'
    (in_fpath, dir_out,
     group_kind, window_width, window_step, value_kind) = set_parameters()
    draw_sequence_distribution(in_fpath, dir_out, group_kind, window_width,
                             window_step, value_kind)

if __name__ == '__main__':
    main()
