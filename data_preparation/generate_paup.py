import re
import glob

files = glob.glob("/Volumes/heloderma4/sonal/encelia/input_files/*nex")

for file in files:
	log = re.sub("nex", "log", file)
	boot = re.sub("nex", "trees", file)
	tre = re.sub("nex", "tre", file)

	out = re.sub("nex", "out", file)
	o = open(out, 'w')
	o.write("Begin paup;\n")
	o.write("set autoclose=yes warntree=no warnreset=no;\n")
	o.write("log start file=%s replace;\n" % log)
	o.write("execute %s;\n" % file)
	o.write("cstatus;\n")
	if re.search("species.nex", file):
		o.write("svdq evalq=all taxpartition=species bootstrap treeFile=%s;\n" % boot)
	else:
		o.write("svdq evalq=all bootstrap treeFile=%s;\n" % boot)
	o.write("savetrees file=%s replace;\n" % tre)
	o.write("end;\n")
	o.write("quit;\n")
	o.close()
	print('/Volumes/heloderma4/sonal/bin/paup4a165_osx %s' % out)
