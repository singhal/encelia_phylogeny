import re
import os
import pandas as pd
import argparse
import random

parser = argparse.ArgumentParser(description="calculate abba-baba")
parser.add_argument('--triad',
                        default=False,
                        help = "triad to test")
parser.add_argument('--out',
                        default=False,
                        help = "outfile")
args = parser.parse_args()
indir = '/Volumes/heloderma4/sonal/encelia/'
snps = os.path.join(indir, 'snps.txt')
inds = os.path.join(indir, 'inds.txt')
sampfile = os.path.join(indir, 'encelia_samples_v4.csv')
d = pd.read_csv(sampfile)
allsp = d.groupby('lineage')['sample'].apply(list).to_dict()

def make_code():
        ix = 0
        code = {}

        allbp = ['A', 'T', 'C', 'G']

        for i, bp1 in enumerate(allbp):
                for bp2 in allbp[i:]:
                        code[str(ix)] = [bp1, bp2]
                        ix += 1
        code['-'] = ['N', 'N']
        return code


def get_inds(inds):
        f = open(inds, 'r')
        allinds = []
        for l in f:
                allinds.append(l.rstrip())
        f.close()
        return allinds


def parse_snp(snps, sp, inds, code, tr, outfile):
	out = tr[-1]
	of = open(outfile, 'w')
	of.write('contig,pos,%s\n' % ','.join(tr))
	
	f = open(snps, 'r')
	for l in f:
		dat = re.split('\s+', l.rstrip())
		genos = dict([(ind, code[s]) for ind, s in zip(inds, dat[2])])
		genos2 = dict([(s, []) for s in sp])
		for s in sp:
			for ind in sp[s]:
				genos2[s] = genos2[s] + genos[ind]
		keep = True

		# check if every species has data
		for s in genos2:
			genos2[s] = [bp for bp in genos2[s] if bp != 'N']
		if min([len(genos2[x]) for x in genos2]) == 0:
			keep = False
 
		# check if outgroup is fixed
		if keep:
			if len(set(genos2[out])) > 1:
				keep = False
		
		# check if two alleles
		g = [genos2[s] for s in genos2]
		g = [bp for g1 in g for bp in g1]
		g = list(set(g))
		if len(g) != 2:
			keep = False

		# okay, this is a winner!
		if keep:
			# ancestral
			a = genos2[out][0]
			# derived
			d = [bp for bp in g if bp != a][0]
			outline = dat[0:2] + [genos2[s].count(d) / float(len(genos2[s])) for s in tr]
			outline = [str(x) for x in outline]
			of.write('%s\n' % ','.join(outline))
	of.close()
	f.close()			
			

code = make_code()
# define the triad
tr = re.split(',', args.triad)
sp = dict([(s, allsp[s]) for s in tr])
inds = get_inds(inds)

# for each line
#	get the sps & inds genos
#	see if it matches qualification
#	calculate af
#	report af
parse_snp(snps, sp, inds, code, tr, args.out)


