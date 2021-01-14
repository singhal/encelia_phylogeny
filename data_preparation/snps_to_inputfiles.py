import re
import operator
import os
import numpy as np
import pandas as pd
import argparse
import random

parser = argparse.ArgumentParser(description="create pop gen files")
parser.add_argument('--structure',
	                default=False,
                	action='store_true',
			help="generate structure input; good for fastStructure and structure")
parser.add_argument('--snapp',
                        default=False,
                        action='store_true',
			help = "generate snapp (nex) file")
parser.add_argument('--svd',
                        default=False,
                        action='store_true',
                        help = "generate svd ind file")
parser.add_argument('--phylonet',
                        default=False,
                        action='store_true',
                        help = "generate phylonet file")
parser.add_argument('--random',
	                default=False,
                	action='store_true',
			help="keep just one snp. default is to use all.")
parser.add_argument('--stem', 
			default=None,
			help="stem for all output files; can just default")
parser.add_argument('--MISS',
			default = None,
			help = "percent missing")
parser.add_argument('--MAC',
			default = None,
			help = "mininum allele count req'd; rec 2")
parser.add_argument('--sample',
                        default = None,
                        help = "# of snps to randomly sample")

args = parser.parse_args()
indir = '/Volumes/heloderma4/sonal/encelia/'
snps = os.path.join(indir, 'snps.txt')
inds = os.path.join(indir, 'inds.txt')
sampfile = os.path.join(indir, 'encelia_samples_v4.csv')
d = pd.read_csv(sampfile)
sp = d.groupby('lineage')['sample'].apply(list).to_dict()

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


def get_af(snp, code):
	a = [code[x] for x in snp]
	a = [x for geno in a for x in geno]
	a = [x for x in a if x != 'N']

	uniq = list(set(a))
	cts = [ a.count(bp) for bp in uniq]

	return len(uniq), min(cts)


def get_random(var):
	var2 = {}
	for c in var:
		keep = random.choice(list(var[c].keys()))
		var2[c] = {}
		var2[c][keep] = var[c][keep]

	return var2


def parse_snps(snps, code, MISS, MAC):
	keep = {}
	f = open(snps, 'r')
	for l in f:
		d = re.split('\s+', l.rstrip())
		snp = d[2]
		miss = snp.count('-')
		if miss / float(len(snp)) < MISS:
			ct, freq = get_af(snp, code)
			if ct == 2 and freq >= MAC:
				if d[0] not in keep:
					keep[d[0]] = {}
				keep[d[0]][d[1]] = snp

	return keep

def get_inds(inds):
	f = open(inds, 'r')
	allinds = []
	for l in f:
		allinds.append(l.rstrip())
	f.close()
	return allinds

def make_snapp(inds, sp, var, stem, code, drop):
	out = stem + '.snapp.nex'
	o = open(out, 'w')

	snp = {	'A': {'A': 'A', 'T': 'W', 'C': 'M', 'G': 'R'},
			'T': {'A': 'W', 'T': 'T', 'C': 'Y', 'G': 'K'},
			'C': {'A': 'M', 'T': 'Y', 'C': 'C', 'G': 'S'},
			'G': {'A': 'R', 'T': 'K', 'C': 'S', 'G': 'G'},
			'N': {'N': '-'}}

	seq = {}
	for ind in inds:
		seq[ind] = ''

	loci = []
	for c in sorted(var.keys()):
		for bp in sorted(var[c].keys()):
			loci.append(c + '.' + str(bp))
			for ind, geno in zip(inds, list(var[c][bp])):
				seq[ind] += snp[code[geno][0]][code[geno][1]]

	for ind in drop:
		del seq[ind]
	# figure out inds to drop
	# only want 1 inds max per species
	keep  = []
	sp2 = {}
	for species, inds in sp.items():
		cts = []
		for ind in inds:
			if ind in seq:
				cts.append((ind, seq[ind].count('-')))
		cts = dict(cts)
		cts = sorted(cts.items(), key=operator.itemgetter(1))
		if len(cts) > 0:
			keep.append(cts[0][0])
			sp2[cts[0][0]] = species
		if len(cts) > 1:
			keep.append(cts[1][0])
			sp2[cts[1][0]] = species

	inds1 = list(seq.keys())
	for ind in inds1:
		if ind not in keep:
			del seq[ind]


	def missing(col):
		return np.count_nonzero(col == '-')

	def mac(col):
		bps = [bp for bp in col if bp != '-']
		uniq = list(set(bps))
		if len(uniq) == 1:
			mac = 0
		else:
			mac = min([bps.count(bp) for bp in uniq])
		return mac

	loci = np.array(loci)
	# because deleted so much data, need to reclean up dataset
	snparr = [list(seq[ind]) for ind in keep]
	snparr = np.array(snparr)
	miss = np.apply_along_axis(missing, 0, snparr)
	miss = miss / float(snparr.shape[0])
	snparr =  snparr[ : , miss < float(args.MISS)]
	loci = loci[ miss < float(args.MISS) ]

	mac = np.apply_along_axis(mac, 0, snparr)
	snparr = snparr[ :, mac >= int(args.MAC)]
	loci = loci[ mac >= int(args.MAC) ] 
	
	c = np.array([re.sub('\..*$', '', locus) for locus in loci])
	uniq_c = list(set(c))
	keeploci = [] 
	for u in uniq_c:
		poss = np.where(c == u)[0]
		keeploci.append(random.sample(poss, 1)[0])
	
	snparr = snparr[ :, keeploci]
	pick = np.random.choice(snparr.shape[1], int(args.sample), replace = False)
	snparr = snparr[ : , pick]

	seq3 = {}
	for ind, snps in zip(keep, snparr):
		seq3[ind] = ''.join(snps)

	maxlen = max([len(ind) for ind in seq]) + 1
	o.write("#NEXUS\n")
	o.write("\n")
	o.write("BEGIN DATA;\n")
	o.write("\tDIMENSIONS NTAX=%s NCHAR=%s;\n" % (len(seq3), len(seq3.values()[0])))
	o.write("\tFORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
	o.write("\tMATRIX\n")
	for ind, s in seq3.items():
		indout = re.sub('-', '_', ind) + '.' + sp2[ind]
		o.write('%s %s\n' % (indout, s))
	o.write(";\nEND;\n")
	o.close()
	

def make_svd(inds, sp, var, stem, code, drop):
	out = stem + '.nex'
	out2 = stem + '.species.nex'
	o = open(out, 'w')
	o2 = open(out2, 'w')

	snp = {	'A': {'A': 'A', 'T': 'W', 'C': 'M', 'G': 'R'},
			'T': {'A': 'W', 'T': 'T', 'C': 'Y', 'G': 'K'},
			'C': {'A': 'M', 'T': 'Y', 'C': 'C', 'G': 'S'},
			'G': {'A': 'R', 'T': 'K', 'C': 'S', 'G': 'G'},
			'N': {'N': '-'}}

	seq = {}
	for ind in inds:
		seq[ind] = ''
	
	for c in sorted(var.keys()):
		for bp in sorted(var[c].keys()):
			for ind, geno in zip(inds, list(var[c][bp])):
				seq[ind] += snp[code[geno][0]][code[geno][1]]

	for ind in drop:
		del seq[ind]
	maxlen = max([len(ind) for ind in seq]) + 1
	
	o.write("#NEXUS\n")
	o2.write("#NEXUS\n")
	o.write("\n")
	o2.write("\n")
	o.write("BEGIN DATA;\n")
	o2.write("BEGIN DATA;\n")
	o.write("\tDIMENSIONS NTAX=%s NCHAR=%s;\n" % (len(seq), len(seq.values()[0])))
	o2.write("\tDIMENSIONS NTAX=%s NCHAR=%s;\n" % (len(seq), len(seq.values()[0])))
	o.write("\tFORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
	o2.write("\tFORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
	o.write("\tMATRIX\n")
	o2.write("\tMATRIX\n")
	for s1 in sorted(sp.keys()):
		for ind in sorted(sp[s1]):
			if ind not in drop:
				s = seq[ind]
				indout = re.sub('-', '_', ind) + ' ' * (maxlen - len(ind))
				o.write('%s %s\n' % (indout, s))
				o2.write('%s %s\n' % (indout, s))
	o.write(";\nEND;\n")
	o2.write(";\nEND;\n")
	o2.write("BEGIN SETS;\n")
	o2.write("\tTAXPARTITION SPECIES =\n")
	ix = 1
	for s1 in sorted(sp.keys()):
		cur = [x for x in sp[s1] if x not in drop]
		if len(cur) > 0:
			o2.write("\t\t%s: %s-%s,\n" % (s1, ix, ix + len(cur) - 1))
			ix = ix + len(cur)
	o2.write(";\nEND;")
	o.close()


def make_structure(inds, b, stem, code, drop):
	o = open(stem + '.str', 'w')

	for ix, ind in enumerate(inds):
		gen1 = [ind, '1', '0', '1', '1']
		gen2 = [ind, '1', '0', '1', '1']

		for c in sorted(b.keys()):
			for pos in sorted(b[c].keys()):
				gen1.append(b[c][pos][ix][0])
				gen2.append(b[c][pos][ix][1])

		if ind not in drop:
			o.write(' '.join(gen1) + '\n')
			o.write(' '.join(gen2) + '\n')

		# if ix == 1:
		#	of.write('num loci: %s\n' % len(gen1) - 5)
		#	of.write('num inds: %s\n' % len(inds) - len(drop))

	o.close()

def make_binary(var, code):
	binary = {}
	for c in var:
		binary[c] = {}
		for pos in var[c]:
			genos = list(var[c][pos])
			genos = [code[geno] for geno in genos]
		
			alleles = [bp for geno in genos for bp in geno]
			alleles = [bp for bp in alleles if bp != 'N']
			alleles = list(set(alleles))
		
			a = {}
			for ix, allele in enumerate(alleles):
				a[allele] = str(ix)
			a['N'] = '-9'

			new = []
			for geno in genos:
				new.append([a[geno[0]], a[geno[1]]])
			binary[c][pos] = new
	return binary

def remove_miss(var, inds):
	tot_len = 0
	indmiss = dict([(ind, 0) for ind in inds])

	for c in var:
		for pos in var[c]:
			tot_len += 1
			for ind, geno in zip(inds, var[c][pos]):
				if geno == '-':
					indmiss[ind] += 1

	to_drop = []
	for ind in inds:
		permiss = indmiss[ind] / float(tot_len)
		if permiss > 0.7:
			to_drop.append(ind)
	
	return to_drop
	
				
def get_sample(var, sample):
	sample = int(sample)
	var2 = {}
	while sample > 0:
		rc = random.choice(var.keys())
		if rc not in var2:
			var2[rc] = {}
		rcpos = random.choice(var[rc].keys())
		if rcpos not in var2[rc]:
			var2[rc][rcpos] = var[rc][rcpos]
			sample = sample - 1
	return var2 

MISS = float(args.MISS)
MAC = int(args.MAC)

code = make_code()
inds = get_inds(inds)
var = parse_snps(snps, code, MISS, MAC)
stem = args.stem + '.miss%s.MAC%s' % (MISS, MAC)

# randomly subsample one SNP per locus
if args.random:
	var = get_random(var)
	stem = stem + '.thinned'
if args.sample and not args.snapp:
	var = get_sample(var, args.sample)
	stem = stem + '.sample%s' % args.sample

drop = remove_miss(var, inds)
out = stem + '.out'
of = open(out, 'w')
for l in drop:
        of.write(l + '\tdrop\n')

# get binary sructure
binary = make_binary(var, code)

if args.structure:
	make_structure(inds, binary, stem, code, drop)

if args.svd:
	make_svd(inds, sp, var, stem, code, drop)

if args.snapp:
	make_snapp(inds, sp, var, stem, code, drop)

of.close()
