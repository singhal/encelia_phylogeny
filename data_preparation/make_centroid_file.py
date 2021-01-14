import re
import os

mfile = '/Volumes/heloderma4/sonal/encelia/clustered/all_encelia.sorted.matches'
seqfile = '/Volumes/heloderma4/sonal/encelia/clustered/all_encelia.sorted.fasta'
outfile = '/Volumes/heloderma4/sonal/encelia/clustered/all_encelia.centroids.fasta'
cutoff = 0.1

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



def write_centroid(m, inds, seq, cutoff):
        for c in m:
                count = len(m[c]) / float(len(inds))
                if count >= cutoff:
			o.write('>%s\n%s\n' % (c, seq[c]))

m, inds = get_matches(mfile)
seq = get_seq(seqfile)
o = open(outfile, 'w')
write_centroid(m, inds, seq, cutoff)
o.close()
