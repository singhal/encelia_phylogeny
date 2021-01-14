import re
import os
import gzip
import glob

indir = '/Volumes/heloderma4/sonal/encelia/ref_variants/'
vcfs = glob.glob(indir + '*gz')

def make_code():
	ix = 0
	code = {}
	for bp1 in ['A', 'T', 'C', 'G']:
		for bp2 in ['A', 'T', 'C', 'G']:
			bps = sorted([bp1, bp2])
			if bps[0] not in code:
				code[bps[0]] = {}
			if bps[1] not in code[bps[0]]:
				code[bps[0]][bps[1]] = ix
				ix = ix + 1
	code['N'] = {}
	code['N']['N'] = '-'
	return code


def get_info(vcfs):
	var = {}

	inds = []
	for vcf in vcfs:
		f = gzip.open(vcf, 'r')
		for l in f:
			if re.search('#CHROM', l):
				d = re.split('\t', l.rstrip())
				inds += d[9:]
				break
		f.close()

	f = gzip.open(vcfs[0], 'r')
	for l in f:
		if re.search('##contig', l):
			c = re.search('##contig=<ID=(\S+),len', l).group(1)
			length = re.search('length=(\d+)', l).group(1)
			var[c] = {}
			for i in range(0, int(length)):
				var[c][str(i + 1)] = '-' * len(inds)
		elif re.search('#CHROM', l):
			break
	f.close()

	return var, inds


def fill_in(var, allinds, code):
	f = gzip.open(vcf, 'r')
	for l in f:
		if re.search('#CHROM', l):
			d = re.split('\t', l.rstrip())
			inds = d[9:]
			inds = [allinds.index(ind) for ind in inds]
		elif not re.search('^#', l) and not re.search('INDEL', l):
                        d = re.split('\t', l.rstrip())
			alleles = [d[3]] + re.split(',', d[4])
			genos = d[9:]

			snpcode = {}
                        for ix, a in enumerate(alleles):
                                snpcode[str(ix)] = a
                        snpcode['.'] = 'N'

			bps = list(var[d[0]][d[1]])
			for ind, geno in zip(inds, genos):
				geno = re.search('(\S\S\S)', geno).group(1)
				geno = re.split('/', geno)
				geno = sorted([snpcode[geno[0]], snpcode[geno[1]]])
				bps[ind] = str(code[geno[0]][geno[1]])
			var[d[0]][d[1]] = ''.join(bps)
	f.close()
	return var

code = make_code()
print(code)
# var, inds = get_info(vcfs)
# for vcf in vcfs:
# 	print(vcf)
# 	var = fill_in(var, inds, code)
# f = open('/Volumes/heloderma4/sonal/encelia/inds.txt', 'w')
# for ind in inds:
# 	f.write('%s\n' % ind)
# f.close()
# f = open('/Volumes/heloderma4/sonal/encelia/snps.txt', 'w')
# for c in var:
#	for pos in var[c]:
# 		f.write('%s %s %s\n' % (c, pos, var[c][pos]))
# f.close()
