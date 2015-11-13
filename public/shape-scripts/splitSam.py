#!/usr/env/bin python
#splitSam.py 11/2015 E.McArthur

# This function splits reads from a sam file and outputs fastq files with reads
# for both the reference and alternate alleles at a particular SNV locus

#USAGE: splitSam.py ref alt chrStr loc random? sam_file_name/list 
# ref = reference allele (A, C, G, T)
# alt = alternate allele (A, C, G, T)
# chrStr = chromosome string - this will be the name of the fasta file you 
#	aligned it to (ex. chr1, ACTB_mRNA)
# loc = location of the SNV
# random? = for reads that don't overlap with the loci: 
#	1 = they will be randomly assigned to the ref and alt fastq files 
#	at the same proportion as the ref and alt alleles
#	0 = reads are thrown away
# sam_file_name/list = either a samfile name or a list of samfile names, the
#	list must be a .txt file with each samfile name on a new row 

#SAMPLE:
# python splitSam.py G A ACTB_mRNA 231 1 M147.sam

#OTHER NOTES: 
# the script will create bam, sorted bam, and bai files for each sam file
# and then deleted when the script finishes running; however, if the file 
# already is created (a bam file with the same name as the sam file) new 
# files won't be created or deleted

import sys
import pysam
import csv
import random
import os
import os.path

refAllele = sys.argv[1]
altAllele = sys.argv[2]
chrStr = sys.argv[3]
loc = int(sys.argv[4]) - 1

if sys.argv[5] == 1:
	random = TRUE

file_name = sys.argv[6]
x_iter = []

if file_name[-4:] == ".txt":
	reader = csv.reader(open(file_name,"r"),delimiter = "\t")
        for rows in reader:
		x_iter.append(rows[0])
else:
	x_iter = [file_name]

for x in x_iter:
	print("File: " + str(x))

	bam =  os.path.isfile(os.path.splitext(x)[0] + ".bam")
	sorted_bam = os.path.isfile(os.path.splitext(x)[0] + ".sorted.bam")
	bai = os.path.isfile(os.path.splitext(x)[0] + ".sorted.bam.bai")

	if not(bam):
		os.system("samtools view -bS " + x + " > " + os.path.splitext(x)[0] + ".bam")
	if not(sorted_bam):
		os.system("samtools sort " + os.path.splitext(x)[0] + ".bam " + os.path.splitext(x)[0] + ".sorted")
	if not(bai):
		os.system("samtools index " + os.path.splitext(x)[0] + ".sorted.bam")

	samfile = pysam.AlignmentFile(os.path.splitext(x)[0] + ".sorted.bam","rb")


	splitRef = []
	splitAlt = []
	count = 0.0
	countRef = 0.0
	countAlt = 0.0
	
	for pileupcolumn in samfile.pileup(chrStr):
		if pileupcolumn.pos == loc:
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del:
					count+=1
					if pileupread.alignment.query_sequence[pileupread.query_position] == refAllele:
						splitRef.append(pileupread.alignment.query_name)
						countRef+=1
					if pileupread.alignment.query_sequence[pileupread.query_position] == altAllele:
                	                        splitAlt.append(pileupread.alignment.query_name)
						countAlt+=1

	
	if not(bam):
		os.system("rm " + os.path.splitext(x)[0] + ".bam")
	if not(sorted_bam):
		os.system("rm " + os.path.splitext(x)[0] + ".sorted.bam")
	if not(bai):
		os.system("rm " + os.path.splitext(x)[0] + ".sorted.bam.bai")
	
	print("\tProportion reference reads = " + str(countRef/count))
	print("\tProportion alternate reads = " + str(countAlt/count))
	print("\tProportion other reads = " + str((count - countRef - countAlt)/count))
	print("\tTotal read count at loci = " + str(count))
	mydict = {}
	reader = csv.reader(open(x,"r"),delimiter = "\t")
	for rows in reader:
		if rows[0][0] != "@":
			if rows[0] not in mydict:
				mydict[rows[0]] = rows[9:11]
			else:
				mydict[rows[0]].extend(rows[9:11])

	f_ref = open(str(os.path.splitext(x)[0]) + "_" + str(chrStr) + "_" + str(loc + 1) + "_ref_" + str(refAllele) + ".fastq","w")
	for i in splitRef:
		if i in mydict:
			tmp = mydict[i]
			while len(tmp) >= 2:
				f_ref.write('@%s\n%s\n+\n%s\n' % (i,tmp[0],tmp[1]))
				tmp = tmp[2:]
			del mydict[i]


	f_alt = open(str(os.path.splitext(x)[0]) + "_" + str(chrStr) + "_" + str(loc + 1) + "_alt_" + str(altAllele) + ".fastq","w")
	for i in splitAlt:
        	if i in mydict:
                	tmp = mydict[i]
                	while len(tmp) >= 2:
                        	f_alt.write('@%s\n%s\n+\n%s\n' % (i,tmp[0],tmp[1]))
                       		tmp = tmp[2:]
                	del mydict[i]
	if random:

		for k,v in mydict.iteritems():
			while len(v) >=2:
				if random.random() >= (countAlt/count):
					f_ref.write('@%s\n%s\n+\n%s\n' % (k,v[0],v[1]))
				else:
					f_alt.write('@%s\n%s\n+\n%s\n' % (k,v[0],v[1]))
				v = v[2:]
	f_ref.close()
	f_alt.close()	
