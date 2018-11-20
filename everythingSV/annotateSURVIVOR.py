
import sys
import re
import os
import subprocess


def split_by_SV_type(avinput, output_prefix, workdir):
		"""
		Splits avinput file by SV type (i.e. in DEL.avinput, DUP.avinput etc)

		Parameters
		----------

		avinput : str
			SVs in avinput format

		output_prefix : str
			sample prefix

		workdir : str
			working directory (output is written here)

   
		Returns
		-------


		Raises
		------

		To do!
		"""
		DEL = workdir + "/" +  output_prefix + ".DEL"
		DEL = open(DEL, 'w')
		DUP = workdir + "/" +  output_prefix + ".DUP"
		DUP = open(DUP, 'w')
		INV = workdir + "/" + output_prefix + ".INV"
		INV = open(INV, 'w')
		INS = workdir + "/" +  output_prefix + ".INS"
		INS = open(INS, 'w')
		BND = workdir + "/" + output_prefix + ".BND"
		BND = open(BND, 'w')
		with open(avinput, 'r') as avinput_file:
			for line in avinput_file:
				line = line.split('\t')
				type = line[5]
				if type == '<DEL>':
					newline="\t".join(line)
					DEL.write("%s\n"%newline)
				elif type == '<DUP>':
					newline="\t".join(line)
					DUP.write("%s\n"%newline)
				elif type == '<INV>':
					newline="\t".join(line)
					INV.write("%s\n"%newline)
				elif type == '<INS>':
					newline="\t".join(line)
					INS.write("%s\n"%newline)
				else:
					newline="\t".join(line)
					BND.write("%s\n"%newline)




def annotate_avinput(avinput, table_annovar, humandb, reference, output_prefix, workdir, *bedfiles):
		"""
		Annotates a set of structural variants in avinput format

		Parameters
		----------

		table_annovar : str
			path to table_annovar.pl script

		humandb : str
			path to humandb/ folder containing annotations of interest

		reference : str
			hg19 or hg38

		output_prefix : str
			sample prefix

		workdir : str
			working directory (output is written here)

		bedfiles : list
			A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, the minimum fraction of overlap of sample SV with region of interest is the third element, and the fourth element is TRUE if it is desired to have the percentage overlap of the database annotation calculated, otherwise FALSE  (e.g. "gene_promoters.bed, 4, 0.9")

		Returns
		-------

		str
			path to annotated SV text file

		Raises
		------

		To do!
		"""


		protocols = "refGene,"
		filenames = ""
		operations = "g,"
		bed_args = ","

		bedfiles = bedfiles[0]
		num_files = len(bedfiles)


		#Iterate through bed files to extract filenames, protocol, and operation for input to Annovar
		for bedfile in bedfiles:
			if num_files > 1:
				bedfile = bedfile.strip('').split(",")
				filenames = filenames + bedfile[0] + ","
				cols = str(bedfile[1])
				frac = str(bedfile[2])
				if "NA" in cols:
					if "NA" not in frac:
						bed_args = bed_args + " -minqueryfrac "  + frac + ","
					else:
						bed_args = bed_args + ","
				elif "NA" not in frac:
					bed_args = bed_args + "-colsWanted " + cols + " -minqueryfrac "  + frac + ","
				else:
					bed_args = bed_args + "-colsWanted " + cols + ","
				protocols = protocols + "bed,"
				operations = operations + "r,"
				num_files = num_files - 1
			else:
				bedfile = bedfile.strip('').split(",")
				filenames = filenames + bedfile[0]
				cols = bedfile[1]
				frac = bedfile[2]
				if "NA" in cols:
					if "NA" not in frac:
						bed_args = bed_args + " -minqueryfrac "  + frac
					else:
						bed_args = bed_args + ""
				elif "NA" not in frac:
					bed_args = bed_args + "-colsWanted " + cols + " -minqueryfrac "  + frac
				else:
					bed_args = bed_args + "-colsWanted " + cols
				protocols = protocols + "bed"
				operations = operations + "r"
				num_files = num_files - 1


		command_line = "perl {} {} {} -buildver {} -out {} -nastring . -remove -otherinfo -protocol {} -bedfile {} -operation {} -arg '{}' ".format(
			table_annovar, avinput, humandb, reference, workdir + '/' +
			output_prefix + '.annovar.header.txt', protocols, filenames, operations, bed_args)

		print command_line

		subprocess.call(command_line, shell=True)

		#Remove annovar header
		command_line = "tail -n +2 {} > {}".format(workdir + '/' +
			output_prefix + '.annovar.header.txt' + '.hg19_multianno.txt', workdir + '/' +
			output_prefix + '.annovar.txt' + '.hg19_multianno.txt')

		subprocess.call(command_line, shell=True)

		return workdir + '/' + output_prefix + ".annovar.txt.hg19_multianno.txt"

def calculate_overlap(output_prefix, workdir, avoutput, svtype, *bedfiles):
		"""
		Takes the *bedfiles argument and extracts the columns of interest where overlap = TRUE, and calculates the fraction of overlap of the sample SV with the database SV/annotation from the  avoutput file

		To do this, the annotation column of the bedfile that is specified must be in the format chr:start:end:<otherinfo>

		Parameters
		----------

		output_prefix : str
			sample prefix

		workdir : str
			working directory (output is written here)

		avoutput: str
			output file from annovar annotation (i.e. <SAMPLE>.annovar.txt.hg19_multianno.txt)

		svtype : str
			DEL, DUP, INV, INS, BND

		bedfiles : list
			A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, the minimum fraction of overlap of sample SV with region of interest is the third element, the fourth element is TRUE if it is desired to have the percentage overlap of the database annotation calculated, otherwise FALSE, and the fifth element is the name of the annotation  (e.g. "gene_promoters.bed, 4, 0.9, TRUE, Gene_promoters")

		Returns
		-------

		str
			path to annotated SV text file

		Raises
		------

		To do!
		"""

		bedfiles = bedfiles[0]

		#The first bed file specified by the user is column 10
		col_num = 10

		#Make a list of column numbers for which a reciprocal overlap calculation is desired
		overlap_col = []
		for bedfile in bedfiles:
			bedfile = bedfile.strip('').split(",")
			overlap = bedfile[3]
			if 'TRUE' in overlap:
				overlap_col.append(col_num)
				col_num = col_num + 1
			else:
				col_num = col_num + 1 

				

		avoutput = open(avoutput, 'r')
		avoutput_overlap_calcs = open(workdir + output_prefix + "." + svtype + ".withoverlap.annovar.hg19_multianno.txt", 'w')

		#Iterate through annovar output file line by line and calculate reciprocal overlap where specified in overlap_col
		for line in avoutput: 
			cols=line.strip('\n').split('\t')
			output_col = []
			#skip header
			if cols[0] == "Chr":
				pass
			else:
				#get sample SV coordinates
				chr=cols[0]
				start=int(cols[1])
				end=int(cols[2])
				length=end-start
				#iterate through each annotation for which an overlap calculation is desired and calculate reciprocal overlap with sample SV
				newcols=[]
				for col in overlap_col:
					annotation = cols[col].split(",")
					recoverlaplist=[]
					for record in annotation:
						if record != '.':
							#AC=record.strip('"').split(":")[3]
							anno_coordinates = record.strip('"').split(":")
							anno_start=int(anno_coordinates[1])
							anno_end=int(anno_coordinates[2])
							if svtype == "INS":
								anno_length = 1
								recoverlap = "1"
							elif svtype == "BND":
								anno_length = 1
								recoverlap = "1"
							else:
								anno_length = float(anno_end-anno_start)
								startMax = max(start, anno_start)
								endMin = min(end, anno_end)
								recoverlap = str((endMin-startMax)/anno_length)
							recoverlaplist.append(recoverlap)
					#create new list associating database annotation coordinates to reciprocal overlap calculation
					new_col =[]
					for anno, rec in zip (annotation, recoverlaplist):
						new_col.append(anno + ":" + rec)
					new_col = ",".join(new_col)
					output_col.append(new_col)
			index = 0
			newcol_index = 0
			final_output_col = []
			#For final output, append columns unchanged where the reciprocal overlap was not calculated. For columns for which a reciprocal overlap was calculated, append modified annotations which include reciprocal overlap fraction
			chr=cols[0]
			start=int(cols[1])
			end=int(cols[2])
			length=str(end-start)
			for col in  cols:
				if cols[0] == "Chr":
					pass
				if index == 3:
					col = length
					final_output_col.append(col)
					index = index + 1
				elif index not in overlap_col: 
					final_output_col.append(col)
					index = index + 1
				else:
					final_output_col.append(output_col[newcol_index])
					newcol_index = newcol_index + 1
					index = index + 1
			newline = newline="\t".join(final_output_col)

			avoutput_overlap_calcs.write("%s\n"%newline)

		return workdir + output_prefix + "." + svtype + ".withoverlap.annovar.hg19_multianno.txt"

					
def add_header(output_prefix, workdir, svtype, *bedfiles):
	 	"""
		Annovar does not label user-provided bed files; this function will add a header with column names to the avoutput file

		Parameters
		----------

		output_prefix : str
			sample prefix

		workdir : str
			working directory (output is written here)

		svtype : str
			DEL, DUP, INV, INS, BND

		bedfiles : list
			A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, the minimum fraction of overlap of sample SV with region of interest is the third element, the fourth element is TRUE if it is desired to have the percentage overlap of the database annotation calculated, otherwise FALSE, and the fifth element is the name of the annotation  (e.g. "gene_promoters.bed, 4, 0.9, TRUE, Gene_promoters")

		-------

		str
			path to annotated SV text file

		Raises
		------

		To do!
		"""   

		header = open(workdir + "/" + output_prefix + "." + svtype +  ".header.txt", 'w')

	 	bedfiles = bedfiles[0]
	 	colnames = ["Chr", "Start", "End", "Length", "Alt", "Func.RefGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene"]
	 	
	 	for bed in bedfiles:
	 		colname = bed.strip('\n').split(',')[4]
	 		colnames.append(colname)

	 	end_colnames = ["SVtype", "Supporting_info", "Format", "LUMPY", "CNVnator", "Manta", "ERDS\n"]
	 	
	 	colnames = colnames + end_colnames

	 	newline = newline="\t".join(colnames)
	 	header.write("%s"%newline)

	 	return workdir + "/" + output_prefix + "." + svtype + ".header.txt"










