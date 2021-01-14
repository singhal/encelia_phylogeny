import re
import subprocess
import os
import argparse
import random
import string
import sys
import pandas as pd

parser = argparse.ArgumentParser(description='Get homology file for groups of interest.')
parser.add_argument('--group', help="group to consider")
parser.add_argument('--c1', help="cluster percentage for ingroup")
parser.add_argument('--c2', help="cluster percentage for outgroup")
parser.add_argument('--out', help="outgroup as inds, comma separated")
args = parser.parse_args()

outdir = '/Volumes/heloderma4/sonal/encelia/homology/'
indfile = '/Volumes/heloderma4/sonal/encelia/encelia_samples_v1.csv'
d = pd.read_csv(indfile)
d = d[d.genome == args.group] 

cldir = '/Volumes/heloderma4/sonal/encelia/cluster_assemblies/'
fadir = '/Volumes/heloderma4/sonal/encelia/ind_assemblies/'

inds = d['sample'].tolist() 
outs = re.split(',', args.out)

orgfiles = [(ind, '%s%s.fasta' % (fadir, ind)) for ind in inds]
orgfiles = dict(orgfiles)
# cluster ingroup seqs at this value
WCLUST1 = float(args.c1)
# cluster outgroup seqs at this value
WCLUST2 = float(args.c2)

inds = ['FAR02']
for ind in inds:
	db = os.path.join(cldir, '%s.fasta' % args.group)
	seq = os.path.join(fadir, '%s.fasta' % ind)
	out1 = os.path.join(outdir, '%s_usearch1' % ind)
	out2 = os.path.join(outdir, '%s_usearch2' % ind)

	if not os.path.isfile(out1):
		subprocess.call("vsearch --usearch_global %s --db %s --id %s --query_cov 0.7 --maxhits 1 --blast6out %s --strand both" % (seq, db, WCLUST1, out1), shell = True)
	if not os.path.isfile(out2):
		subprocess.call("vsearch --usearch_global %s --db %s --id %s --maxhits 1 --blast6out %s --strand both" % (db, seq, WCLUST1, out2), shell = True)

	matches1 = {}
	matches2 = {}
	m1 = {}
	m2 = {}

	f = open(out1, 'r')
	for l in f:
		d = re.split('\t', l.rstrip())
		strand = '+'
		if int(d[7]) < int(d[6]):
			strand = '-'
		percent = float(d[2])

		matches1[d[0]] = {'match': d[1], 'strand': strand, 'per': percent, 'recip': 'NA'}
		if d[1] not in m1:
			m1[d[1]] = []
		m1[d[1]].append(d[0])
	f.close()

	f = open(out2, 'r')
        for l in f:
		d = re.split('\t', l.rstrip())
		matches2[d[0]] = d[1]
	
		if d[1] not in m2:
                        m2[d[1]] = []
                m2[d[1]].append(d[0])
	f.close()
	
	for seq in matches1:
		match = matches1[seq]['match']

		if seq in m2 and match in m1:
			if len(m2[seq]) == 1 and len(m1[match]) == 1:
				if matches2[match] == seq:
					matches1[seq]['recip'] = 'EASY'
				else:
					matches1[seq]['recip'] = 'DITCH'
			elif len(m2[seq]) == 1 and len(m1[match]) > 1:
				pass
			elif len(m2[seq]) > 1 and len(m1[match]) == 1:
				pass
			else:
				pass
		else:
			if len(m1[match]) == 1:
				matches1[seq]['recip'] = 'KEEP'
				print(seq)
			else:
				pass

