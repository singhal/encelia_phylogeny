import re
import glob
import os
import argparse

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('-c','--cov', help='cov', required=True)
args = parser.parse_args()
cov = args.cov

all = {}

def get_seq(seqfile):
        id = ''
        seq = {}
        
        f = open(seqfile, 'r')
        for l in f:
                if re.search('>', l):
                        id = re.search('>(\S+)', l).group(1)
                        seq[id] = ''
                else:
                        seq[id] += l.rstrip()

        f.close()
	return seq

files = glob.glob('/Volumes/heloderma4/sonal/encelia/variable_PRG/*cov_%s*fasta' % cov)
for ix, file in enumerate(files):
	print(ix)
	ind = re.search('([^/]+)\.cov_', file).group(1)
	seq = get_seq(file)

	for c, s in seq.items():
		if c not in all:
			all[c] = {}
		all[c][ind] = s

outdir = '/Volumes/heloderma4/sonal/encelia/phylogeny_v2/alignments_cov%s/' % cov
for c in all:
	out = os.path.join(outdir, '%s.fasta.aln' % c)
	o = open(out, 'w')
	for id, s in all[c].items():
		o.write('>%s\n%s\n' % (id, s))
	o.close()
