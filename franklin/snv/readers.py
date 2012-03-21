'''
Created on 2011 aza 21

@author: peio
'''
class VcfParser(object):
    'A vcf reader'
    def __init__(self, fpath):
        'Class initiator'
        self._fpath = fpath
        self.header = None

        if float(self.version) < 4:
            self._get_header_v3()
        else:
            self._get_header_v4()
        self._index = None

    def _get_version(self):
        'version of the vcf'
        fhand = open(self._fpath)
        first_line = fhand.readline()
        fhand.close()
        first_line = first_line.strip()
        first_line =  first_line.lstrip('##')
        kind, value = first_line.split('=', 1)
        if kind not in ('fileformat', 'format'):
            msg = 'Vcf file malformed. first line must be version info'
            raise ValueError(msg)
        return value.split('v')[1]

    version = property(_get_version)

    def _get_header_v3(self):
        'it returns the header'
        if self.header is not None:
            return self.header
        headers = {}
        for line in open(self._fpath):
            if not line.startswith('#'):
                break
            if line.startswith('##'):
                line = line.strip()
                line = line.lstrip('##')
                kind, value = line.split('=', 1)
                if kind == 'FILTER':
                    if kind not in headers:
                        headers[kind] = {}
                    filter_type, filter_info = value.split(',', 1)
                    headers[kind][filter_type] = filter_info.strip('"')
                elif kind in ('FORMAT', 'INFO'):
                    if kind not in headers:
                        headers[kind] = {}
                    name, example, type_, desc = value.split(',')
                    headers[kind][name] = {'type':type_, 'example':example,
                                           'description':desc}
                else:
                    headers[kind] = value
            else:
                line = line.lstrip('#')
                headers['colnames'] = line.split()
        self.header = headers

    def _get_header_v4(self):
        'it returns the header'
        if self.header is not None:
            return self.header
        headers = {}
        for line in open(self._fpath):
            if not line.startswith('#'):
                break
            if line.startswith('##'):
                line = line.strip()
                line = line.lstrip('##')
                kind, value = line.split('=', 1)
                if kind == 'FILTER':
                    if kind not in headers:
                        headers[kind] = {}
                    values = self._parse_header_line_v4(value)
                    headers[kind][values['ID']] = values['Description']
                elif kind in ('FORMAT', 'INFO'):
                    if kind not in headers:
                        headers[kind] = {}
                    values = self._parse_header_line_v4(value)
                    name = values['ID']
                    del values['ID']
                    headers[kind][name] = values
            else:
                line = line.lstrip('#')
                headers['colnames'] = line.split()

        self.header = headers

    @staticmethod
    def _parse_header_line_v4(line):
        'it parses a header line, but olny the value part. between <>'
        values = {}
        line = line.strip('<')
        line = line.strip('>')
        for pair in line.split(','):
            key, value = pair.split('=')
            values[key] = value.strip('"')
        return values

    def _get_vcfs(self):
        'vcf generator'
        colnames = self.header['colnames']
        for line in open(self._fpath):
            if line.startswith('#'):
                continue
            yield self._parse_vcf_line(line, colnames)
    vcfs = property(_get_vcfs)

    def _parse_vcf_line(self, line, colnames):
        '''It parses the cvf svn line'''
        vcf_items = line.split()
        vcf = dict(zip(colnames[:9], vcf_items[:9]))
        # reformat FILTER
        vcf['FILTER'] = vcf['FILTER'].split(';')
        # REformat INFO
        info = vcf['INFO']
        vcf['INFO'] = {}
        for info_ in info.split(';'):
            info_key, info_value = info_.split('=')
            vcf['INFO'][info_key] = info_value
        # reformat FORMAT
        format_string = vcf['FORMAT']
        vcf['FORMAT'] = {}
        for format_ in format_string.split(';'):
            format_key, format_value = format_.split(':')
            vcf['FORMAT'][format_key] = format_value


        vcf['samples'] = {}
        for samples in zip(colnames[9:], vcf_items[9:]):
            allele_count = {}
            alleles, values =  samples[1].split(':')
            for index , allele in enumerate(alleles.split('|')):
                if allele in (',', '.'):
                    continue
                if allele == '0':
                    allele = vcf['REF']
                else:
                    allele = vcf['ALT'].split(',')[int(allele) - 1]
                try:
                    count_ = int(values.split(',')[index])
                except ValueError:
                    continue
                allele_count[allele] = count_
            vcf['samples'][samples[0]] = allele_count
        return vcf

    def _make_index(self):
        '''it makes an index of the vcf file. It takes the vcf position
        (chrom, position) as index'''
        if self._index is not None:
            return self._index
        index = {}
        fhand = open(self._fpath, 'rt')
        rawline = 'filled'
        while len(rawline) != 0:
            prior_tell = fhand.tell()
            rawline = fhand.readline()
            if  rawline and rawline[0] == '#':
                continue
            index[tuple(rawline.split()[:2])] = prior_tell
        self._index = index

    def get_snv(self, position):
        'It returns an snv giving it position'
        colnames = self.header['colnames']
        if self._index is None:
            self._make_index()
        fhand = open(self._fpath)
        file_position = self._index[position]
        fhand.seek(file_position)
        return self._parse_vcf_line(fhand.readline(), colnames)



