import argparse
import glob
import os
import re
import subprocess

"""
Sonal Singhal
created on 23 June 2016
Written assuming nothing!
"""

def get_args():
        parser = argparse.ArgumentParser(
                        description="This creates the files that then get " 
                                    "aligned in the next script.",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter
                        )
   
        # dir
        parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help='directory with alignments'
                )

        # output dir
        parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for concatenated file'
                )

	# outfile stem
	parser.add_argument(
                '--stem',
                type=str,
                default=None,
                help='stem for output files'
                )

	# inds to drop
	parser.add_argument(
		'--drop',
		type=str,
		default=None,
		help='inds to drop'
		)

	return parser.parse_args()

def get_seq(seq):
	f = open(seq, 'r')
	id = ''
	s = {}
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			# get rid of reverse
			id = re.sub('^_R_', '', id)
			# get rid of ending
			id = re.sub('_\d+$', '', id)			

			s[id] = ''
		else:
			s[id] += l.rstrip()

	f.close()

	loc_length = len(s[s.keys()[0]])
	name = re.search('([^/]+).fasta', seq).group(1)

	return s, loc_length, name


def make_concatenated(args):
	outdir = args.outdir

	out1 = os.path.join(outdir, '%s.concatenated.phy' % args.stem)
	out2 = os.path.join(outdir, '%s.concatenated.partitions' % args.stem)

	files = glob.glob(args.dir + '/*aln')

	all = {}
	loci = {}
	inds = {}

	for locus in files:
		seq, loc_length, name = get_seq(locus)

		loci[name] = loc_length
		all[name] = seq			

		for ind in seq:
			if ind not in inds:
				inds[ind] = 1

	if args.drop:
		drop = re.split(',', args.drop)
	else:
		drop = []
	inds = [ind for ind in inds if ind not in drop]

	f1 = open(out1, 'w')
	f2 = open(out2, 'w')
	totseq = sum(loci.values())
	f1.write('%s\t%s\n' % (len(inds), totseq))
	for ind in inds:
		s = ''
		for locus in sorted(loci.keys()):
			if ind in all[locus]:
				s += all[locus][ind]
			else:
				s += '-' * loci[locus]
		f1.write('%s\t%s\n' % (ind, s))
	f1.close()

	cur_pos = 1
	for locus in sorted(loci.keys()):
		end_pos = cur_pos + loci[locus] - 1
		f2.write('%s\t%s\t%s\n' % (locus, cur_pos, end_pos))
		cur_pos = end_pos + 1
	f2.close()
		

def main():
	args = get_args()
	make_concatenated(args)

if __name__ == "__main__":
	main()
