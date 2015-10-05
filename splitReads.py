#!/usr/env/bin python
#USAGE: splitReads.py A C 0.48 chr1
#USAGE: splitReads.py ref alt mutRate chrStr

import sys
import pysam
import csv
import random

refAllele = sys.argv[1]
altAllele = sys.argv[2]
mutRate = float(sys.argv[3])
chrStr = sys.argv[4]
loc = int(sys.argv[5]) - 1


x_iter = ["PR","MR","DC"]

for x in x_iter:
	samfile = pysam.AlignmentFile("/home/emcarthur/transcriptome/out" + x + ".sorted.bam","rb")


	splitRef = []
	splitAlt = []

	for pileupcolumn in samfile.pileup(chrStr):
		if pileupcolumn.pos == loc:
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del:
					if pileupread.alignment.query_sequence[pileupread.query_position] == refAllele:
						splitRef.append(pileupread.alignment.query_name)
				if pileupread.alignment.query_sequence[pileupread.query_position] == altAllele:
                	                        splitAlt.append(pileupread.alignment.query_name)


	mydict = {}
	reader = csv.reader(open("/home/emcarthur/transcriptome/out" + x + ".sam","r"),delimiter = "\t")
	for rows in reader:
		if len(rows) > 11:
			if rows[0] not in mydict:
				mydict[rows[0]] = rows[9:11]
			else:
				mydict[rows[0]].extend(rows[9:11])
		 
	f_ref = open("/home/emcarthur/transcriptome/" + x + "ref.fastq","w")
	for i in splitRef:
		if i in mydict:
			tmp = mydict[i]
			while len(tmp) >= 2:
				f_ref.write('@%s\n%s\n+\n%s\n' % (i,tmp[0],tmp[1]))
				tmp = tmp[2:]
			del mydict[i]


	f_alt = open("/home/emcarthur/transcriptome/" + x + "alt.fastq","w")
	for i in splitAlt:
        	if i in mydict:
                	tmp = mydict[i]
                	while len(tmp) >= 2:
                        	f_alt.write('@%s\n%s\n+\n%s\n' % (i,tmp[0],tmp[1]))
                       		tmp = tmp[2:]
                	del mydict[i]

	for k,v in mydict.iteritems():
		while len(v) >=2:
			if random.random() >= mutRate:
				f_ref.write('@%s\n%s\n+\n%s\n' % (k,v[0],v[1]))
			else:
				f_alt.write('@%s\n%s\n+\n%s\n' % (k,v[0],v[1]))
			v = v[2:]
	f_ref.close()
	f_alt.close()	
