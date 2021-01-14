import argparse
import glob
import os
import re
import subprocess
import gzip

"""
Sonal Singhal
created on 31 May 2019
Written assuming:
	* samtools 1.3.1
	* picard 2.4.1
	* bwa
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Call SNPs by individual. Written assuming "
                            " samtools 1.3.1, picard 2.4.1, and"
                            " bwa",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# sample
	parser.add_argument(
                '--sample',
                type=str,
                default=None,
                help='sample for which to run script.'
                )

        # bwa
	parser.add_argument(
                '--bwa',
                type=str,
                default=None,
                help='bwa executable, full path.'
                )

	# samtools
        parser.add_argument(
                '--samtools',
                type=str,
                default=None,
                help='samtools executable, full path.'
                )

        # bcftools
        parser.add_argument(
                '--bcftools',
                type=str,
                default=None,
                help='bcftools executable, full path.'
                )


	# picard
        parser.add_argument(
                '--picard',
                type=str,
                default=None,
                help='picard executable, full path.'
                )
	
	# CPUs
	parser.add_argument(
                '--CPU',
                type=int,
                default=1,
                help='# of CPUs to use in alignment.'
               )

	# memory
        parser.add_argument(
                '--mem',
                type=int,
                default=1,
                help='Memory available, as an int, in terms of Gb.'
               )
               
        # outdir
	parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for alignments.'
                )

	# readdir
	parser.add_argument(
                '--readdir',
                type=str,
                default=None,
                help="Full path to files with reads."
                )

	# PRG
	parser.add_argument(
                '--prg',
                type=str,
                default=None,
                help="Full path to pseudoref genome."
                )

	return parser.parse_args()


def get_info(args):
	# get the genome
	genome = args.prg
	s = args.sample
	outdir = args.outdir
    
	# makes the outdir
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	reads = {}
	reads[s] = [os.path.join(args.readdir, '%s_R1.final.fq.gz' % s),
			os.path.join(args.readdir, '%s_R1.final.fq.gz' % s),
			os.path.join(args.readdir, '%s_unpaired.final.fq.gz' % s)]

	return reads, genome, outdir


def prepare_seq(args, genome):
	# does all the prep necessary for the PRG
	if not os.path.isfile(genome + '.fai'):
		subprocess.call("%s faidx %s" % (args.samtools, genome), shell=True)
	out = re.sub('.fa.*', '.dict', genome)
	if not os.path.isfile(out):
		subprocess.call("%s CreateSequenceDictionary R=%s O=%s" % 
                                (args.picard, genome, out), shell=True)
	if not os.path.isfile(genome + '.bwt'):
		subprocess.call("%s index %s" % (args.bwa, genome), shell = True)


def align_seq(args, sample, r, seq):
	dir = args.outdir

	out1 = os.path.join(dir, '%s.sam' % sample)
	out1b = os.path.join(dir, '%s_u.sam' % sample)
	out2 = os.path.join(dir, '%s.mateFixed.bam' % sample)
	out3a = os.path.join(dir, '%s_p.mateFixed.sorted.bam' %  sample)
	out3b = os.path.join(dir, '%s_u.sorted.bam' %  sample)
	out3 = os.path.join(dir, '%s.mateFixed.sorted.bam' %  sample)
	out4 = os.path.join(dir, '%s.rg.mateFixed.sorted.bam' % (sample))

	# need a tmpdir for when sorting BAM files
	tmpdir = os.path.join(dir, sample)
	if not os.path.isdir(tmpdir):
		os.mkdir(tmpdir)

	# align using bwa
	subprocess.call("%s mem -t %s %s %s %s > %s" % (args.bwa, args.CPU, seq, r[0], r[1], out1), shell=True)
	subprocess.call("%s mem -t %s %s %s > %s" % (args.bwa, args.CPU, seq, r[2], out1b), shell=True)
	# fixmate
        subprocess.call("%s FixMateInformation I=%s O=%s" % (args.picard, out1, out2), shell=True)
        # sorted
        subprocess.call("%s sort -O bam -o %s -T %s %s" % (args.samtools, out3a, tmpdir, out2), shell=True)
        subprocess.call("%s sort -O bam -o %s -T %s %s" % (args.samtools, out3b, tmpdir, out1b), shell=True)
        # merge
        subprocess.call("%s merge %s %s %s" % (args.samtools, out3, out3a, out3b), shell=True)
        # readgroup
        subprocess.call("%s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % (args.picard, out3, out4, sample, sample, sample), shell=True)
        subprocess.call("%s index %s" % (args.samtools, out4), shell=True)

	# not marking duplicates because RADseq expects many duplicates
	# not doing indel realigner because too slow
	
	# remove the files
	[os.remove(x) for x in [out1, out1b, out2, out3, out3a, out3b]]
	
	# remove the dir
	os.rmdir(tmpdir)

	return out4


def call_snps(args, files, genome):
	outdir = args.outdir

        out1 = os.path.join(outdir, '%s.bcf' % (args.sample))
        out2 = os.path.join(outdir, '%s.vcf.gz' % (args.sample))

	if not os.path.isfile(out2):
        	subprocess.call("%s mpileup -ABIug -t AD -f %s -o %s %s" % (args.samtools, genome, out1, ' '.join(files)), shell=True)
		# output variant sites only
		subprocess.call("%s call -mO z -o %s %s" % (args.bcftools, out2, out1), shell=True)

	os.remove(out1)
    
def main():
	# get arguments
	args = get_args()
	reads, genome, outdir = get_info(args)

	inds = sorted(list(reads.keys()))
	# round 1
	# prep sequence
	prepare_seq(args, genome)
	# do the alignments
	bamfiles1 = []
	for ind in inds:
		bamout = align_seq(args, ind, reads[ind], genome)
		bamfiles1.append(bamout)
	call_snps(args, bamfiles1, genome)


if __name__ == "__main__":
	main()
