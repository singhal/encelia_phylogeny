import argparse
import collections
from ete3 import Tree
import gzip
import itertools as it
import numpy as np
import os
import pandas as pd
import pickle
import random
import re

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tree", required=True, help='Newick string.')
parser.add_argument("-o", "--outgroup", required=True, help='Outgroup.')
args = parser.parse_args()


def get_triads(sps, tree, args):
	# first need to identify all legal three groups
	in_sps = [x for x in sps if x != args.outgroup]
	triads = []
	for group in it.combinations(in_sps, 3):
		cur_min = 1e6
		sp1 = None
		sp2 = None

		# and sort them
		for ix, x in enumerate(group):
			for y in group[(ix+1):]:
				dist = tree.get_distance(x, y)
				if dist < cur_min:
					sp1 = x
					sp2 = y
					cur_min = dist
		left = [x for x in group if x not in [sp1, sp2]]
		# and the constant outgroup
		triad = [sp1, sp2, left[0], args.outgroup]
		triads.append(triad)
	return triads

tree = Tree(args.tree)
sps = [leaf.name for leaf in tree]
triads = get_triads(sps, tree, args)
for triad in triads:
	print(','.join(triad))