#!/usr/env/bin python
# This function splits reads from a sam file at multiple loci and outputs fastq files with reads
# for both the reference and alternate alleles all specified SNV loci

#USAGE: splitSam_multLoc.py mutFile.txt random? sam_file_name/list
# mutFile.txt = a .txt tab-separated file that has each mutaion Loci in its own row
#	the first column should the be reference allele, the second column should be
# 	the alternate allele, the third column should be the chromosome string (ex. chr1,
# 	ACTB_mRNA), and the fourth column should be the numerical location of the SNV
# 	that you want to split on
#	A sample row would be: "G\tA\tACTB_mRNA\t231"
# random? = for reads that don't overlap with the loci: 
#	1 = they will be randomly assigned to the ref and alt fastq files 
#	at the same proportion as the ref and alt alleles
#	0 = reads are thrown away
# sam_file_name/list = either a samfile name or a list of samfile names, the
#	list must be a .txt file with each samfile name on a new row 

#SAMPLE:
# python mutFile.txt 1 M147.sam

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

mutLists = []
reader = csv.reader(open(sys.argv[1]),delimiter="\t")
for rows in reader:
	mutLists.append(rows)

if sys.argv[2] == 1:
	random = TRUE

file_name = sys.argv[3]
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

	for mutList in mutLists:
		splitRef = []
		splitAlt = []
		count = 0.0
		countRef = 0.0
		countAlt = 0.0
		
		for pileupcolumn in samfile.pileup(mutList[2]):
			if pileupcolumn.pos == int(mutList[3]) - 1:
				for pileupread in pileupcolumn.pileups:
					if not pileupread.is_del:
						count+=1
						if pileupread.alignment.query_sequence[pileupread.query_position] == mutList[0]:
							splitRef.append(pileupread.alignment.query_name)
							countRef+=1
						if pileupread.alignment.query_sequence[pileupread.query_position] == mutList[1]:
       	         	                       		splitAlt.append(pileupread.alignment.query_name)
							countAlt+=1
	
		print("Loci: " + str(mutList[3]))
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
	
		f_ref = open(str(os.path.splitext(x)[0]) + "_" + str(mutList[2]) + "_" + str(mutList[3]) + "_ref_" + str(mutList[0]) + ".fastq","w")
		for i in splitRef:
			if i in mydict:
				tmp = mydict[i]
				while len(tmp) >= 2:
					f_ref.write('@%s\n%s\n+\n%s\n' % (i,tmp[0],tmp[1]))
					tmp = tmp[2:]
				del mydict[i]
	
	
		f_alt = open(str(os.path.splitext(x)[0]) + "_" + str(mutList[2]) + "_" + str(mutList[3]) + "_alt_" + str(mutList[1]) + ".fastq","w")
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

	if not(bam):
                os.system("rm " + os.path.splitext(x)[0] + ".bam")
        if not(sorted_bam):
                os.system("rm " + os.path.splitext(x)[0] + ".sorted.bam")
        if not(bai):
                os.system("rm " + os.path.splitext(x)[0] + ".sorted.bam.bai")	
