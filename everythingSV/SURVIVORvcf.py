
import sys
import re
import os
import subprocess

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

	def SURVIVOR_to_avinput(self, output_prefix, workdir, size):
		"""
		Takes the vcf output of SURVIVOR and converts into a form compatible with annovar (avinput) i.e. chr start end ref alt otherinfo, where ref and alt are represented by "0" for compatibility

		Parameters
		----------

		output_prefix : str
			prefix for output file

		workdir : str
			working directory (output is written here)

		size: str
			'small' (50-400000bp) or 'large' (5000+bp)

		Returns
		-------

		str
			path to avinput file

		Raises
		------

		To do!

		"""
		SURVIVORvcf = open(self.file_path, 'r')
		avinput = open(workdir + "/" + output_prefix + ".SURVIVOR." + size + ".avinput", 'w')

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

		return workdir + "/" + output_prefix + ".SURVIVOR." + size + ".avinput"





