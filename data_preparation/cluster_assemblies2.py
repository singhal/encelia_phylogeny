import re
import os

mfile = '/Volumes/heloderma4/sonal/encelia/clustered/all_encelia.sorted.matches'
alndir = '/Volumes/heloderma4/sonal/encelia/phylogeny/alignments/'
seqfile = '/Volumes/heloderma4/sonal/encelia/clustered/all_encelia.sorted.fasta'
cutoff = 0.2

def get_matches(mfile):
	m = {}
	inds = {}

	f = open(mfile, 'r')
	for l in f:
		d = re.split('\t', l.rstrip())
		ind1 = re.sub('_\d+$', '', d[0])
		ind2 = re.sub('_\d+$', '', d[1])
		if ind1 not in inds:
			inds[ind1] = 1
		if ind2 not in inds:
			inds[ind2] = 1

		if ind1 != ind2:
			if d[1] not in m:
				m[d[1]] = {}
			m[d[1]][d[0]] = d[2]
	f.close()

	return m, inds


def get_seq(seqfile):
	id = ''
        seq = {}
        f = open(seqfile, 'r')

        for l in f:
                if re.search('>', l):
                        id = re.search('>(\S+)', l.rstrip()).group(1)
                        seq[id] = ''
                else:
                        seq[id] += l.rstrip().upper()

        f.close()

        return seq


def rev_comp(dna):
	complement = {'N': 'N', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return ''.join([complement[base] for base in dna[::-1]])


def write_alns(m, inds, seq, cutoff):
	for c in m:
		count = len(m[c]) / float(len(inds))
		if count >= cutoff:
			o = open(os.path.join(alndir, '%s.fasta' % c), 'w')
			o.write('>%s\n%s\n' % (c, seq[c]))
			for c1, orr in m[c].items():
				s = seq[c1]
				if orr == '-':
					s = rev_comp(s)
				o.write('>%s\n%s\n' % (c1, s))
			o.close()


m, inds = get_matches(mfile)
seq = get_seq(seqfile)
write_alns(m, inds, seq, cutoff)
