import re
import gzip
import argparse

# get sample name
def get_args():
        parser = argparse.ArgumentParser(
                description="Create pseudo-reference genome",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
                )

	# sample
        parser.add_argument(
                '--sample',
                type=str,
                default=None,
                help='sample for which to run script.'
                )

	# min cov
	parser.add_argument(
		'--mincov',
		type = int,
		default = None,
		help = 'minimum coverage required to keep variant'
		)

	return parser.parse_args()


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

        for id, s in seq.items():
                seq[id] = len(s) * ['-']

        return seq


def get_var(seq, var, cov):
	f = gzip.open(var, 'r')
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l.rstrip())
			# check quality
			if float(d[5]) >= 20:
				# check coverage
				bpcov = int(re.search('DP=(\d+)', d[7]).group(1))
				if bpcov >= cov:
					pos = int(d[1]) - 1
					c = d[0]
				
					ref = d[3]
					if re.search('^0/0', d[9]):
						seq[c][pos] = ref
					else:
						alleles = [d[3]] + re.split(',', d[4])
						genos = re.search('^(\S\S\S)', d[9]).group(1)
						genos = [alleles[int(x)] for x in re.split('/', genos)]
						seq[c][pos] = codes[genos[0]][genos[1]]
	f.close()
	return seq						

codes = {'A': {'A': 'A', 'T': 'W', 'C': 'M', 'G': 'R'}, 
        'T': {'A': 'W', 'T': 'T', 'C': 'Y', 'G': 'K'},
        'C': {'A': 'M', 'T': 'Y', 'C': 'C', 'G': 'S'},
        'G': {'A': 'R', 'T': 'K', 'C': 'S', 'G': 'G'}}

args = get_args()
# get ref genome
ref = '/Volumes/heloderma4/sonal/encelia/ref_bias/%s.fasta' % args.sample
seq = get_seq(ref)

# get variants
var = '/Volumes/heloderma4/sonal/encelia/variants/%s.vcf.gz' % args.sample
seq = get_var(seq, var, int(args.mincov))

# do output
out = '/Volumes/heloderma4/sonal/encelia/variable_PRG/%s.cov_%s.fasta' % (args.sample, args.mincov)
o = open(out, 'w')
for c, bp in seq.items():
	miss = bp.count('-') / float(len(bp))
	if miss < 0.5:
		o.write('>%s\n%s\n' % (c, ''.join(bp)))
o.close()
	
