import re
import glob
import pandas as pd
import os
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description="Run final alignment steps.")
parser.add_argument('--cl', help="Cluster for which to run script.")
parser.add_argument('--file', help="file to use for sample data.")
parser.add_argument('--dir', help="out directory.")
args = parser.parse_args()
cl = args.cl
c_file1 = args.file

min_qual = 20

d = pd.read_csv(c_file1)
d = d.groupby('lineage')
clusters = dict([(cluster, sorted(group['sample'].tolist())) for cluster, group in d])
inds = clusters[cl]
print(len(inds))

gatk = '/Volumes/heloderma4/sonal/bin/GenomeAnalysisTK.jar'
samtools = '/Volumes/heloderma4/sonal/bin/samtools-1.3.1/samtools'
bcftools = '/Volumes/heloderma4/sonal/bin/bcftools/bcftools'

dir = os.path.join(args.dir, 'ref_alignments')
seqdir = os.path.join(args.dir, 'ref_bias') 

bamfiles = [os.path.join(dir, '%s.rg.mateFixed.sorted.bam' %  ind) for ind in inds]
bamfiles= [x for x in bamfiles if os.path.isfile(x)]	

raw_vcf = os.path.join(dir, '%s.raw.vcf' % cl)
filt_vcf = os.path.join(dir, '%s.filtered.vcf' % cl)
seq = os.path.join(seqdir, '%s.fasta' % cl)

subprocess.call("%s mpileup -ugf %s %s | %s call -vmO v -o %s" % (samtools, seq, ' '.join(bamfiles), bcftools, raw_vcf), shell=True)

# filter the vcf
f = open(raw_vcf, 'r')
o = open(filt_vcf, 'w')

for l in f:
	if re.match('^#', l):
		o.write(l)
	else:
		d = re.split('\t', l.rstrip())
		if float(d[5]) >= min_qual:
			o.write(l)
f.close()
o.close()
