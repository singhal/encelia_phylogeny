import argparse
import glob
import os
import re
import subprocess
import gzip

"""
Sonal Singhal
created on 22 March 2017
Written assuming:
	* samtools 1.3.1
	* picard 2.4.1
	* bwa
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Reference correct genome. Written assuming "
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


def align_seq(args, sample, r, seq, ix):
	dir = args.outdir

	out1 = os.path.join(dir, '%s.sam' % sample)
	out1b = os.path.join(dir, '%s_u.sam' % sample)
	out2 = os.path.join(dir, '%s.mateFixed.bam' % sample)
	out3a = os.path.join(dir, '%s_p.mateFixed.sorted.bam' %  sample)
	out3b = os.path.join(dir, '%s_u.sorted.bam' %  sample)
	out3 = os.path.join(dir, '%s.mateFixed.sorted.bam' %  sample)
	out4 = os.path.join(dir, '%s.rg.mateFixed.sorted%s.bam' % (sample,ix))

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
                seq[id] = list(s)

        return seq


def call_snps(args, files, genome, ix):
	outdir = args.outdir

        out1 = os.path.join(outdir, '%s.%s.bcf' % (args.sample, ix))
        out2 = os.path.join(outdir, '%s.%s.vcf.gz' % (args.sample, ix))

	if not os.path.isfile(out2):
        	subprocess.call("%s mpileup -ABIug -f %s -o %s %s" % (args.samtools, genome, out1, ' '.join(files)), shell=True)
		# output variant sites only
		subprocess.call("%s call -mvO z -o %s %s" % (args.bcftools, out2, out1), shell=True)

        # get in seq
        seq = get_seq(genome)

        # mutate seq
        i = gzip.open(out2, 'r')
        for l in i:
                if not re.search('#', l):
                        d = re.split('\t', l.rstrip())
                        # check biallelic
                        if d[3] in ['A', 'T', 'C', 'G'] and d[4] in ['A', 'T', 'C', 'G']:
                                # get af
                                genos = d[9:]
				genos = [re.search('(\S/\S)', x).group(1) for x in genos]
                                genos = [re.split('/', x) for x in genos]
                                genos = [x for ind in genos for x in ind]
                                genos = [x for x in genos if x != '.']

                                if len(genos) > 0:
                                        af = genos.count('1') / float(len(genos))
                                        # mutate seq
                                        if af >= 0.5:
                                                pos = int(d[1]) - 1
                                                # print(seq[d[0]][pos])
						seq[d[0]][pos] = d[4]
						# print(seq[d[0]][pos])
						# print('***')
						
        i.close()

        out = os.path.join(outdir, '%s.%s.fasta' % (args.sample, ix))
        o = open(out, 'w')
        for id, s in seq.items():
                o.write('>%s\n%s\n' % (id, ''.join(s)))
        o.close()


	os.remove(out1)
	os.remove(out2)
        return out


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
		bamout = align_seq(args, ind, reads[ind], genome, 1)
		bamfiles1.append(bamout)
	genome1 = call_snps(args, bamfiles1, genome, 1)

        # round2
        # prep sequence
        prepare_seq(args, genome1)
        # do the alignments for the second time
        bamfiles2 = []
        for ind in inds:
                bamout = align_seq(args, ind, reads[ind], genome1, 2)
                bamfiles2.append(bamout)
        # call the variants for the second time
        genome2 = call_snps(args, bamfiles2, genome1, 2)
        
        # round3
        # prep sequence
        prepare_seq(args, genome2)
        # do the alignments for the third time
        bamfiles3 = []
        for ind in inds:
                bamout = align_seq(args, ind, reads[ind], genome2, 3)
                bamfiles3.append(bamout)
        # call the variants for the third time
        genome3 = call_snps(args, bamfiles3, genome2, 3)

	# round4
        # prep sequence
        prepare_seq(args, genome3)
        # do the alignments for the fourth time
        bamfiles4 = []
        for ind in inds:
                bamout = align_seq(args, ind, reads[ind], genome3, 4)
                bamfiles4.append(bamout)
        # call the variants for the third time
        genome4 = call_snps(args, bamfiles4, genome3, 4)
	final = re.sub('\.4\.', '.', genome4)
	subprocess.call("cp %s %s" % (genome4, final), shell=True)

	files1 = glob.glob(args.outdir + '%s.1.*' % args.sample)
	files2 = glob.glob(args.outdir + '%s.2.*' % args.sample)
	files3 = glob.glob(args.outdir + '%s.3.*' % args.sample)
	files4 = glob.glob(args.outdir + '%s.4.*' % args.sample)
	files5 = glob.glob(args.outdir + '%s*bam*' % args.sample) 
	rmfiles = files1 + files2 + files3 + files4 + files5
	
	[os.remove(x) for x in rmfiles]
	

if __name__ == "__main__":
	main()
