#Step 1: SURVIVOR vcf to avinput (custom python script)
#Step 2: Separate by SV type
#Step 3: add 1 bp to INS
#Step 4: annovar for each SV type
#Step 5: remove header, concatenate (just remove header from all but first)
#Step 6: reciprocal overlap calculations: DGV, 1000G, DECIPHER, ISCA
#Step 7: remove #name
#Step 8: reciprocal overlap with inhouse (CAUSES and IMAGINE
#Step 9: add header
#Step 10: remove 'score'
#Step 11: filter out SUPP=1, inhouse count etc



#vcf module
#SURVIVOR to SV2
#SURVIVOR to avinput, separate by type, add 1 bp to INS
#annovar
#overlap calculations for everything, remove NAME and SCORE


#to do: get a table annotator working
import sys
import re
import os
import subprocess
import pandas as pd

class SURVIVOR_vcf():

	"""

	A class for handling and annotating SURVIVOR files

	Attributes
	----------

	file_path: str
		path to SURVIVOR vcf

	"""

	def __init__(self, file_path):
		self.file_path = file_path

	def SURVIVOR_to_avinput(self, output_prefix, workdir):
		"""
		Takes the vcf output of SURVIVOR and converts into a form compatible with annovar (avinput) i.e. chr start end ref alt otherinfo, where ref and alt are represented by "0" for compatibility

		Parameters
		----------

		output_prefix : str
			prefix for output file

		workdir : str
			working directory (output is written here)

		Returns
		-------

		str
			path to avinput file

		Raises
		------

		To do!

		"""
		SURVIVORvcf = open(self.file_path, 'r')
		avinput = open(workdir + "/" + output_prefix + ".SURVIVOR.avinput", 'w')

		for line in SURVIVORvcf:
			if line.startswith('#'):
				pass
			else:
				cols=line.strip('\n').split('\t')
				chr = cols[0]
				ref = "0"
				alt = "0"
				sv = cols[4]
				other = cols[7]
				geno = cols [8:]
				if sv == '<INV>':
					start = other.split(';')
					start = start[6].split('=')
					start = int (start[1])
					end = int(str(cols[1]))
					if start < end:
						start = start
					else:
						end = start
						start = str(cols[1])
				elif ":" in sv:
					start = str(cols[1])
					end = int(cols[1]) + 1
					end = str(end)
				else:
					start = str(cols[1])
					end = other.split(';')
					end = end[6].split('=')
					end = end[1]
				colsnew = []
				colsnew.extend((chr,str(start),str(end),ref,alt,sv,other))
				colsnew = colsnew + cols[8:]
				newline="\t".join(colsnew)
				avinput.write("%s\n"%newline)

		return workdir + "/" + output_prefix + ".SURVIVOR.avinput"





	def annotate_avinput(self, table_annovar, humandb, reference, output_prefix, workdir, *bedfiles):
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

		bed_file : str
			regions of interest in bed format

		cols_wanted : int
			column number for annotation of interest

		minqueryfrac : float
			minimum fraction of overlap of sample SV with region of interest for annotation

		output_prefix : str
			sample prefix

		workdir : str
			working directory (output is written here)

		bedfiles : list
			A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, and the minimum fraction of overlap of sample SV with region of interest is the third element (e.g. "gene_promoters.bed, 4, 0.9")

		cols_wanted : int
			column number for annotation of interest

		minqueryfrac : float
			minimum fraction of overlap of sample SV with region of interest for annotation

		Returns
		-------

		str
			path to annotated SV text file

		Raises
		------

		To do!

		"""
		protocols = ""
		filenames = ""
		operations = ""
		bed_args = ""
		num_files = len(bedfiles)

		#Iterate through bed files to extract filenames, protocol, and operation for input to Annovar
		for bedfile in bedfiles:
			if num_files > 1:
				bedfile = bedfile.strip('').split(",")
				filenames = filenames + bedfile[0] + ","
				cols = str(bedfile[1])
				frac = str(bedfile[2])
				print cols
				print frac
				if "NA" in cols:
					print cols
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
					print cols
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
			table_annovar, self.file_path, humandb, reference, workdir + '/' +
			output_prefix + '.annovar.txt', protocols, filenames, operations, bed_args)

		print command_line

		subprocess.call(command_line, shell=True)

		return command_line

#SURVIVORvcf = SURVIVOR_vcf("/home/mcouse/projects/rrg-frid/IMAGINE/SV/BATCH4/CP038-P.SURVIVOR.SV.200.vcf")

#avinput = SURVIVORvcf.SURVIVOR_to_avinput("CP038-P","/home/mcouse/projects/rrg-frid/IMAGINE/SOFTWARE/StructuralVariantAnalysis/everythingSV/")

#avinput = SURVIVOR_vcf("/home/mcouse/projects/rrg-frid/IMAGINE/SV/BATCH4/CP038-P.SURVIVOR.SV.200.avinput")

#avinput.annotate_avinput ("/home/mcouse/projects/rrg-frid/IMAGINE/SOFTWARE/annovar/table_annovar.pl",
#{}"/home/mcouse/projects/rrg-frid/IMAGINE/SOFTWARE/annovar//humandb",
#{}"hg19", "CP038-P", #"/home/mcouse/projects/rrg-frid/IMAGINE/SOFTWARE/StructuralVariantAnalysis/ever#ythingSV/", "InHouseDB_SV_50_20180507_DUP_name.bed, NA, NA", #"InHouseDB_SV_50_20180507_DEL_name.bed, NA, NA", "InHouseDB_SV_50_20180507_INV_name.bed, 8, NA",
#{}"InHouseDB_SV_50_20180507_INS_name.bed, 8, 0.8" )
