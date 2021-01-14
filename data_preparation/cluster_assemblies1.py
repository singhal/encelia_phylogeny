import re
import subprocess
import glob
import os
import argparse

parser = argparse.ArgumentParser(description='reads')
parser.add_argument('--ind', help="ind to run this on")
args = parser.parse_args()

WCLUST = 0.9
vsearch = 'vsearch'

ind = os.path.join('/Volumes/heloderma4/sonal/encelia/ind_assemblies2/', '%s.fasta' % args.ind)

indtmp1 = ind + '_1'
indtmp2 = ind + '_2'
indtmp3 = re.sub('.fasta', '.clustered.fasta', ind)
if not os.path.isfile(indtmp3):
	subprocess.call("%s --derep_fulllength %s --output %s --fasta_width 0 --strand both" % (vsearch, ind, indtmp1), shell=True)
	subprocess.call("%s --sortbylength %s --output %s" % (vsearch, indtmp1, indtmp2), shell=True)
	subprocess.call("%s --cluster_smallmem %s --centroids %s --id %s --usersort --fasta_width 0 --strand both --minsl 0.5 --query_cov 0.7" % (vsearch, indtmp2, indtmp3, WCLUST), shell=True)
	os.remove(indtmp1)
	os.remove(indtmp2)
