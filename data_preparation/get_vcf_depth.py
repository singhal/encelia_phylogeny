import re
import os
import subprocess
import argparse
import pandas as pd
import gzip

parser = argparse.ArgumentParser(description='Get raw VCF and depth file for clusters.')
parser.add_argument('--cl', help="cluster to run this on")
parser.add_argument('--dir', help="directory to run this on")
args = parser.parse_args()
cl = args.cl
maindir = args.dir

c_file = os.path.join(maindir, 'encelia_samples_v4.csv')
seq_dir = os.path.join(maindir, 'ref_bias')
bam_dir = os.path.join(maindir, 'ref_alignments')
vcf_dir = os.path.join(maindir, 'ref_variants')
cov_dir = vcf_dir
samtools = '/Volumes/heloderma4/sonal/bin/samtools-1.3.1/samtools'
bcftools = '/Volumes/heloderma4/sonal/bin/bcftools/bcftools'

min_qual = 20
min_cov = 3
# define cap as anything 3 times over the median (xs_cov)
xs_cov = 10

##################
# prep output dirs
##################
if not os.path.isdir(vcf_dir):
	os.mkdir(vcf_dir)
if not os.path.isdir(cov_dir):
	os.mkdir(cov_dir)

def get_clusters(c_file):
        d = pd.read_csv(c_file)
        d = d.groupby('lineage')

	clusters = dict([(name, sorted(g['sample'].tolist()))  for name, g in d])

        return clusters

def get_raw_vcf(vcf_dir, bam_dir, seq_dir, cl, inds):
	raw_bcf = os.path.join(vcf_dir, '%s.bcf' % cl)
	raw_vcf = os.path.join(vcf_dir, '%s.raw.vcf.gz' % cl)
	bam_files = [os.path.join(bam_dir, '%s.rg.mateFixed.sorted.bam') % ind for ind in inds]
	seq = os.path.join(seq_dir, '%s.fasta' % cl)	

	subprocess.call("%s mpileup -ABIg -t AD -f %s -o %s %s" % (samtools, seq, raw_bcf, ' '.join(bam_files)), shell=True)
	subprocess.call("%s call -mO z -o %s %s" % (bcftools, raw_vcf, raw_bcf), shell=True)

	return raw_vcf
	

def get_depth(cov_dir, vcf, cl, inds):
	depthout = os.path.join(cov_dir, '%s.depth_profiling.csv' % cl)
	o = open(depthout, 'w')
	o.write('cluster,individual,depth,counts\n')

	f = gzip.open(vcf, 'r')	

	depth = {}
	for ind in inds:
		depth[ind] = {}
	
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l.rstrip())
			
			tags = re.split(':', d[8])
			# identify where AD is
			ad = tags.index('AD')

			if float(d[5]) >= min_qual:
				genos = d[9:]
				for ind, geno in zip(inds, genos):
					geno = re.split(':', geno)
					ind_ad = re.split(',', geno[ad])
					ind_ad = sum([int(x) for x in ind_ad])

					if ind_ad > 0:
						if ind_ad not in depth[ind]:
							depth[ind][ind_ad] = 0
						depth[ind][ind_ad] += 1
	f.close()

	for ind in depth:
		for d, count in depth[ind].items():
			o.write('%s,%s,%s,%s\n' % (cl, ind, d, count))
	o.close()

	return depthout

def summarize_depth(depthout):
	d = pd.read_csv(depthout)
	g = d.groupby('individual')

	maxcov = {}
	for name, d in g:
		d = d.sort_values(by='depth')
		d['depth_sum'] = d.counts.cumsum()
		midpoint = int(d.depth_sum.max() / 2.0)
		median = d[d.depth_sum >= midpoint].depth.tolist()[0]

                cap = median * xs_cov
		if cap < 20:
			cap = 20

                maxcov[name] = cap
		print(name, cap)
	return maxcov

def filter_vcf(cl, inds, raw_vcf, vcf_dir, maxcov):
	filt_vcf = os.path.join(vcf_dir, '%s.qual_filtered%s.cov_filtered_min%s_max%sx.vcf.gz' % (cl, min_qual, min_cov, xs_cov))

	f = gzip.open(raw_vcf, 'r')
	o = gzip.open(filt_vcf, 'w')

	for l in f:
		if re.match('^#', l):
			o.write(l)
		else:
			d = re.split('\t', l.rstrip())

			tags = re.split(':', d[8])
			# identify where AD is
			ad = tags.index('AD')

			if float(d[5]) >= min_qual:
				genos = d[9:]
				num_miss = 0
				for ix, (id, geno) in enumerate(zip(inds, genos)):
					ind = re.split(':', geno)
                                        ind_ad = re.split(',', ind[ad])
                                        ind_ad = sum([int(x) for x in ind_ad])
                                        # depth for ind is below what we allow
                                        if ind_ad < min_cov or ind_ad > maxcov[id]:
                                                num_miss += 1
                                                d[9+ix] = re.sub('^\S\/\S', './.', d[9+ix]) 
                                # don't print genotype if missing in all inds
                                if num_miss < len(genos):
                                        o.write('\t'.join(d) + '\n')
                                        
        f.close()
        o.close()


clusters = get_clusters(c_file)
raw_vcf = get_raw_vcf(vcf_dir, bam_dir, seq_dir, cl, clusters[cl])
depthout = get_depth(cov_dir, raw_vcf, cl, clusters[cl])
maxcov = summarize_depth(depthout)
filter_vcf(cl, clusters[cl], raw_vcf, vcf_dir, maxcov)
