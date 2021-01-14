import re 
import glob
import pandas as pd
import os
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description="Run final alignment steps.")
parser.add_argument('--ind', help="individual for which to run sript")
parser.add_argument('--file', help="file to use for sample data.")
parser.add_argument('--dir', help="out directory.")
args = parser.parse_args()

ind = args.ind
c_file1 = args.file

d = pd.read_csv(c_file1)
cl = d.ix[d['sample'] == ind, 'lineage'].tolist()[0]

gatk = '/Volumes/heloderma4/sonal/bin/GenomeAnalysisTK.jar'
samtools = '/Volumes/heloderma4/sonal/bin/samtools-1.3.1/samtools'
bcftools = '/Volumes/heloderma4/sonal/bin/bcftools/bcftools'

dir = os.path.join(args.dir, 'ref_alignments')
seq_dir = os.path.join(args.dir, 'ref_bias')

orig = os.path.join(dir, '%s.rg.mateFixed.sorted.bam' %  ind)
filt_vcf = os.path.join(dir, '%s.filtered.vcf' % cl)
seq = os.path.join(seq_dir, '%s.fasta' % cl)

final = os.path.join(dir, '%s.rg.mateFixed.sorted.final.bam' % ind)
if not os.path.isfile(final):
	recal = os.path.join(dir, '%s.recal.table' % ind)

	subprocess.call("%s index %s" % (samtools, orig), shell=True)
	subprocess.call('java -Xmx10g -jar %s -T BaseRecalibrator -R %s -knownSites %s -I %s -o %s -nct 4' % (gatk, seq, filt_vcf, orig, recal), shell=True)
	subprocess.call('java -Xmx10g -jar %s -T PrintReads -R %s -I %s --BQSR %s -o %s ' % (gatk, seq, orig, recal, final), shell=True) 

	os.remove(recal)
	# os.remove(orig)
