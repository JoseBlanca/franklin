'Blast utilities'
from biolib.biolib_utils import call, create_temp_fasta_files

def _parse_tabular_bl2seq(output):
    'It parses a tabular bl2seq output'
    hsps = []
    for line in output.splitlines():
        hsp = {}
        if line.startswith('#'):
            continue
        items = line.split()
        name1 = items[0]
        name2 = items[1]
        hsp['evalue'] = float(items[10])
        hsp['score'] = int(items[11])
        hsp['alignment'] = {}
        hsp['alignment'][name1] = {}
        hsp['alignment'][name2] = {}
        hsp['alignment'][name1]['start'] = int(items[6]) - 1
        hsp['alignment'][name1]['end'] = int(items[7]) - 1
        hsp['alignment'][name2]['start'] = int(items[8]) - 1
        hsp['alignment'][name2]['end'] = int(items[9]) - 1 
        hsps.append(hsp)
    return hsps

def bl2seq(seq1, seq2, evalue=1e-10, program='blastn'):
    'It does a bl2seq an it returns the result.'
    fileh1, fileh2 = create_temp_fasta_files(seq1, seq2)
    filen1 = fileh1.name
    filen2 = fileh2.name
    #we run the blast
    cmd = ['bl2seq', '-i', filen1, '-j', filen2, '-p', program,
           '-e', str(evalue), '-m', 'T', '-D', '1']
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('Problem running bl2seq: '+ stderr)
    fileh1.close()
    fileh2.close()
    result = _parse_tabular_bl2seq(stdout)
    return result

def _parse_water(output):
    'It parses the water output'
    ali_section = False
    score = None
    n_line_in_ali = 0
    ali_lines = [None, None, None, None]
    for line in output.splitlines():
        if 'Score:' in line:
            score = float(line.split()[-1])
        if score and (not line or line.isspace()):
            ali_section = True
            continue
        if ali_section:
            #bad line
            if line.startswith('#'):
                continue
            #bad line
            try:
                if not line.split()[1].isalnum():
                    continue
            except:
                continue
            if n_line_in_ali == 0:
                ali_lines[0] = line #fist line of the alignment
            if n_line_in_ali == 1:
                ali_lines[1] = line #second line of the alignment
            ali_lines[2] = ali_lines[3] #line before last of the alignment
            ali_lines[3] = line  #last line of the alignment
            n_line_in_ali += 1
    #now we get the start and end
    result = {}
    result['score'] = score
    result['alignment'] = {}
    name1 = ali_lines[0].split()[0]
    result['alignment'][name1] = {}
    name2 = ali_lines[1].split()[0]
    result['alignment'][name2] = {}
    result['alignment'][name1]['start'] = int(ali_lines[0].split()[1]) - 1
    result['alignment'][name2]['start'] = int(ali_lines[1].split()[1]) - 1
    result['alignment'][name1]['end'] = int(ali_lines[2].split()[-1]) - 1
    result['alignment'][name2]['end'] = int(ali_lines[3].split()[-1]) - 1
    return result

def water(seq1, seq2, gapopen=20):
    'It does a water alignment an it returns the result.'
    fileh1, fileh2 = create_temp_fasta_files(seq1, seq2)
    filen1 = fileh1.name
    filen2 = fileh2.name
    #we run the blast
    cmd = ['water', filen1, filen2, '-stdout', '-auto', '-snucleotide1',
           '-snucleotide2', '-gapopen', str(gapopen)]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('Problem running water: '+ stderr)
    fileh1.close()
    fileh2.close()
    result = _parse_water(stdout)
    return result 
