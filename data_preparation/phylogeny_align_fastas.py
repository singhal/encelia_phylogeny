import glob
import os

alns = glob.glob('/Volumes/heloderma4/sonal/encelia/phylogeny/alignments/*fasta')
for aln in alns:
	out = aln + '.aln'
	if not os.path.isfile(out):
		print('mafft %s > %s.aln' % (aln, aln))
