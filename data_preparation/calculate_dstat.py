import argparse
import collections
import itertools as it
import pickle
import gzip
import numpy as np
import os
import random
import re

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--triad", required=True, help='triad string')
parser.add_argument("-n", "--name", required=True, help='Output name')
args = parser.parse_args()

indir = '/Volumes/heloderma4/sonal/encelia/'
outfile = os.path.join(indir, 'pop_gen', 'dstat_%s.csv' % args.name)
n_boot = 200

def get_variants(indir, sps):
	var = {}

	for sp in sps:
		print(sp)
		vcf = os.path.join(indir, 'ref_variants', 
			'%s.qual_filtered20.cov_filtered_min3_max10x.vcf.gz' % sp)
		f = gzip.open(vcf, 'r')

		for l in f:
			if not re.search('#', l) and not re.search('INDEL', l):
				d = re.split('\t', l.rstrip())

				c = d[0]
				pos = int(d[1])
				alleles = [d[3]] + re.split(',', d[4])

				if c not in var:
					var[c] = {}
				if pos not in var[c]:
					var[c][pos] = {}
				var[c][pos][sp] = []
				
				snpcode = {}
				for ix, a in enumerate(alleles):
					snpcode[str(ix)] = a
				snpcode['.'] = 'N'

				genos = d[9:]
				genos = [re.search('^(\S/\S)', x).group(1) \
							for x in genos]
				for geno in genos:
					geno = re.split('/', geno)
					snps = [snpcode[g] for g in geno]
					var[c][pos][sp] += snps

	return var


def get_sp_allele(var2):
	alleles = var2
	alleles = [a for a in alleles if a != 'N']

	allele = None
	# needs to be rep'd by at least 2 chrs in each sp
	if len(alleles) < 2:
		allele = 'N'
	# no polymorphic sites
	elif len(set(alleles)) > 1:
		allele = 'N'
	else:
		allele = alleles[0]

	return(allele)

def get_snps(triad, var):
	d_stat = {}

	for exon in var:
		d_stat[exon] = {'abba': 0, 'baba': 0}
		for pos in var[exon]:
			# get alleles
			aln = []
			for sp in triad:
				if sp in var[exon][pos]:
					aln.append(get_sp_allele(var[exon][pos][sp]))
				else:
					aln.append('N')

			keep = True
			# don't keep any with undefined
			if aln.count('N') == 0:
				counts = collections.Counter(aln)
				for a in counts:
					# only biallelics 
					# with each one rep'd twice
					if counts[a] != 2:
						keep = False
			else:
				keep = False
					
			if keep:
				if aln[0] != aln[1]:
					if aln[0] == aln[2] and aln[1] == aln[3]:
						d_stat[exon]['baba'] += 1
					elif aln[0] == aln[3] and aln[1] == aln[2]:
						d_stat[exon]['abba'] += 1
	
	return d_stat
				
def calc_dstat(dstat):
	exons = list(dstat.keys())
	top = np.sum([dstat[c]['abba'] - dstat[c]['baba'] for c in exons])
	bottom = np.sum([dstat[c]['abba'] + dstat[c]['baba'] for c in exons])

	# number of informative markers
	num_inf = np.sum([1 for c in exons if dstat[c]['abba'] > 0 or dstat[c]['baba'] > 0])

	# dstat
	if bottom != 0:
		dval = top / float(bottom)
	else:
		dval = np.nan

	abba = np.sum([dstat[c]['abba'] for c in exons])
	baba = np.sum([dstat[c]['baba'] for c in exons])

	return dval, num_inf, abba, baba

def sample_wr(population, k):         
    "used for bootstrap sampling"
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack
    ix = [_int(_random() * n) for i in it.repeat(None, k)]
    which = [population[i] for i in ix]

    return which

def calc_boot(dstat, times):
	exons = list(dstat.keys())
	boot = []
	for ix in range(times):
		sample = sample_wr(exons, len(exons))

		top = np.sum([dstat[c]['abba'] - dstat[c]['baba'] for c in sample])
		bottom = np.sum([dstat[c]['abba'] + dstat[c]['baba'] for c in sample])

		# dstat
		if bottom != 0:
			dval = top / float(bottom)
		else:
			dval = np.nan

		boot.append(dval)

	return boot


out = open(outfile, 'w')
out.write('sp1,sp2,sp3,outgroup,dstat,boot_std,zscore,num_inf_loci,abba,baba\n')
sps = re.split(',', args.triad)
# get variant data
# what should my data structure be?
#exon, pos, inds, variants?
# var = get_variants(indir, sps)
# pickle.dump(var, open( "%s_variants.pickle" % args.name, "wb" ))
var = pickle.load(open( "%s_variants.pickle" % args.name, "rb" ))
# print(sps)
# print(inds)
dstat = get_snps(sps, var)
full_dstat, num_inf_loci, abba, baba = calc_dstat(dstat)
boot = calc_boot(dstat, n_boot)
boot_std = np.std(boot)
Z = abs(full_dstat / float(boot_std))
out.write('%s,%s,%s,%s,%s,%s,%s\n' % (','.join(sps), full_dstat, boot_std, Z, num_inf_loci, abba, baba))
out.close()
