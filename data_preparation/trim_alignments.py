import re
import glob
import os
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--seqlen", required=True, help='minimum seq len to retain.')
parser.add_argument("-n", "--perloci", required=True, help='minimum percent loci for ind to retain.')
parser.add_argument("-i", "--perinds", required=True, help='minimum percent inds for loci to retain.')
parser.add_argument("-m", "--miss", required=True, help='minimum completeness to retain.')
parser.add_argument("-d", "--dir", required=True, help='alignment dir.')
parser.add_argument("-o", "--out", required=True, help='outfile stem')
args = parser.parse_args()

seqlendrop = int(args.seqlen)
underrepresented = float(args.perloci)
inds_needed = float(args.perinds)
missing = float(args.miss)
aln_dir = args.dir
mainout = args.out

outdir = re.sub('alignments', 'trim_%s_%s_%s_%s' % (seqlendrop, underrepresented, inds_needed, missing), aln_dir)
if not os.path.isdir(outdir):
	os.mkdir(outdir)
mainout = mainout + '_%s_%s_%s_%s.out' % (seqlendrop, underrepresented, inds_needed, missing)
mo = open(mainout, 'w')

def get_seq(x):
	f = open(x, 'r')
	id = ''
	seq = {}
	for l in f:
		if re.search('^>(\S+)', l.rstrip()):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			id = re.sub('^_R_', '', id)
			seq[id] = ''
		else:
			seq[id] += l.rstrip().upper()

	return seq

alns = glob.glob(aln_dir + '/*aln')
print(len(alns))

badinds = []
badloci = {}
inds = {}

for ix, aln in enumerate(alns):
	locname = re.search('([^/]+)\.fasta', aln).group(1)
	seqaln = get_seq(aln)
	for ind in seqaln:
		ind = re.sub('_\d+$', '', ind)
		if ind not in inds:
			inds[ind] = 0
		inds[ind] += 1
	if ix % 100 == 0:
	 	print(ix)

for bad in badinds:
        mo.write('%s,ALL,USER_DELETED\n' % bad)

# inds: too underrepresented
num_loci = len(alns)
for ind in inds:
	percent = inds[ind] / float(num_loci)
	if percent < underrepresented:
		badinds.append(ind)
		mo.write('%s,ALL,TOO_UNDERREPRESENTED\n' % (ind))

inds = [ind for ind in inds if ind not in badinds] 

# columns with too much missing
for ix, aln in enumerate(alns):
        locname = re.search('([^/]+)\.fasta', aln).group(1)
        seqaln = get_seq(aln)

	ids = list(seqaln.keys())
	for badind in badinds:
		bad = [id for id in ids if re.search('%s_\d+$' % badind, id)]	
		if len(bad) > 0:
			del seqaln[bad[0]]

	for ind, s in seqaln.items():
		slen = len(s) - s.count('-') - s.count('N') - s.count('?')
		if slen < seqlendrop:
			del seqaln[ind]
			mo.write('%s,%s,TOO_SHORT\n' % (ind, locname))

	if len(seqaln) / float(len(inds)) >= inds_needed:
		seq = [True] * len(list(seqaln.values())[0])
		# for each column
		for ix, pos in enumerate(seq):
			# identify all sites
			sites = [s[ix] for c, s in seqaln.items()]
			# identify portion missing
			permiss = len(sites) - (sites.count('-') + sites.count('N') + sites.count('?'))
			permiss = permiss / float(len(sites))
			# identify if it is good or not
			if permiss < missing:
				seq[ix] = False
		out = os.path.join(outdir, '%s.fasta.aln' % locname)
		o = open(out, 'w')
		for c, s in seqaln.items():
			newseq = ''.join([bp if keep else '' for bp, keep in zip(s, seq)])
			newseq = newseq.upper()
			o.write('>%s\n%s\n' % (c, newseq))
		o.close()
		mo.write('%s,%s,NUM_SITES_REMOVED\n' % (seq.count(False), locname))

	else:
		mo.write("%s,ALL,TOO_FEW_INDS\n" % (locname))
