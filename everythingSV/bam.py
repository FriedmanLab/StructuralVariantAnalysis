import sys
import os
import subprocess
import pandas as pd

"""
Note: need to run this script in python v 2.7, not 3 or greater
"""
class bam_file():

	"""

	A class for handling and processing bam files

	Attributes
	----------

	file_path: str
		path to bam file

	reference: str
		path to reference fasta

	samtools_path: str
		path to samtools executable

	"""

	def __init__(self, file_path, reference, samtools_path):
		self.file_path = file_path
		self.reference = reference
		self.samtools_path = samtools_path

	def extract_discordants(self, output_prefix):
		"""
		Runs LUMPY on a bwa-mem aligned whole genome bam to detect structural
		variants.

		Parameters
		----------

		reference : str
			path to reference genome
		lumpy_dir : str
			path to LUMPY directory
		samtools_path : str
			path to samtools executable
		output_prefix: str
			path and prefix for output files

		Returns
		-------

		str
			path to unsorted discordant reads bam

		Raises
		------

		To do!

		"""
		command_line = "samtools view -b -F 1294 {} > {}.discordants.unsorted.bam".format(
				self.file_path, output_prefix)
		subprocess.call(command_line, shell=True)
		return output_prefix


	def extract_splitters(self, lumpy_dir, output_prefix):
		"""
		Extracts split reads from bwa-mem aligned whole genome bam.

		Parameters
		----------

		lumpy_dir : str
				path to LUMPY directory
		samtools_path : str
				path to samtools executable
		working_dir : str
				path to working directory
		output_prefix: str
				path and prefix for output files

		Returns
		-------
		str
			path to unsorted split reads bam

		Raises
		------

		To do!

		"""
		command_line = "samtools view -h  {} | {}/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - >  {}.splitters.unsorted.bam".format(
				self.file_path, lumpy_dir, output_prefix)
		subprocess.call(command_line, shell=True)
		return output_prefix

	def sort_bam(self, bam_path, output_prefix):
		"""
		Sort bam file

		Parameters
		----------

		bam_path: str
				path to bam file

		output_prefix: str
				path and prefix for output files

		Returns
		-------
		str
			path to sorted bam file

		Raises
		------

		To do!

		"""
		command_line = "samtools sort {} -o {}.sorted.bam -@ 2".format(
				bam_path, output_prefix)
		subprocess.call(command_line, shell=True)
		return output_prefix



	def get_insert_size(self, lumpy_dir, output_prefix):
		"""
		Get read insert size

		Parameters
		----------

		lumpy_dir : str
				path to lumpy directory

		output_prefix: str
				path and prefix for output files

		insert_size: str
				path and prefix for insert size file

		Returns
		-------
		list
			[mean insert size, standard deviation of insert size]

		Raises
		------

		To do!

		"""
		command_line = "samtools view {} | tail -n+100000 | {}/scripts/pairend_distro.py -r 150 -X 4 -N 10000 -o {} > {}  ".format(
				self.file_path, lumpy_dir, output_prefix + ".histo", output_prefix +  ".insertsize")
		subprocess.call(command_line, shell=True)
		insert_stats = output_prefix + ".insertsize"
		insert_stats = pd.read_csv(insert_stats, sep='\t', names=["Mean", "SD"])
		mean = float(insert_stats.iloc[0,0].split(":")[1])
		sd = float(insert_stats.iloc[0,1].split(":")[1])
		insert_stats = [mean, sd]


		return insert_stats


	def run_LUMPY(
		self, lumpy, output_prefix, insert_stats, read_length, min_mapq, num_support):
		"""
		Run LUMPY on bwa-mem aligned bam file to detect SVs

		Parameters
		----------

		lumpy : str
				path to lumpy executable

		output_prefix: str
				path and prefix for output files

		insert_stats: list
				output of get_insert_size()

		read_length: int
				length of read in bp

		min_mapq: int
				minimum mapping quality of reads to be used

		num_support: int
				minimum number of reads that support an SV


		Returns
		-------
		str
			path to LUMPY vcf output

		Raises
		------

		To do!

		"""
		ins_mean = insert_stats[0]
		ins_sd = insert_stats[1]

		command_line = "{} -mw {} -tt 0 -pe id:{},bam_file:{},histo_file:{},mean:{},stdev:{},read_length:{},min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,min_mapping_threshold:{} -sr id:{},bam_file:{},back_distance:20,weight:1,min_mapping_threshold:{} > {}".format(
			lumpy, num_support, output_prefix, output_prefix + ".discordants.sorted.bam", output_prefix + ".insertsize.histo",
		ins_mean, ins_sd, read_length,min_mapq,output_prefix,output_prefix + ".splitters.sorted.bam",
		min_mapq,output_prefix + ".LUMPY.vcf")

		subprocess.call(command_line, shell=True)
		return output_prefix + ".LUMPY.vcf"

	def run_CNVnator(
		self, chromosomes, output_prefix):
		"""
		Runs CNVnator on a bwa-mem aligned whole genome bam to detect CNVs.

		Parameters
		----------

		chromosomes : str
			path to directory with chromosome.fa files
		output_prefix: str
			path and prefix for output files

		Returns
		-------

		str
			path to CNVnator vcf

		Raises
		------

		To do!

		"""
		command_line = "cnvnator -root {} -genome {} -tree {}".format(
			output_prefix + ".root", self.reference, self.file_path)
		subprocess.call(command_line, shell=True)
		command_line = "cnvnator -root {} -genome {} -his {} -d {}".format(
			output_prefix + ".root", self.reference, 1000, chromosomes)
		subprocess.call(command_line, shell=True)
		command_line = "cnvnator -root {} -genome {} -stat {}".format(
			output_prefix + ".root", self.reference, 1000)
		subprocess.call(command_line, shell=True)
		command_line = "cnvnator -root {} -genome {} -partition {}".format(
			output_prefix + ".root", self.reference, 1000)
		subprocess.call(command_line, shell=True)
		command_line = "cnvnator -root {} -genome {} -call {} > {}".format(
			output_prefix + ".root", self.reference, 1000, output_prefix + ".CNVcall.1000")
		subprocess.call(command_line, shell=True)
		command_line = "cnvnator2VCF.pl {} {} > {}".format (
			output_prefix + ".CNVcall.1000", chromosomes, output_prefix + ".CNVcall.1000.vcf" )
		subprocess.call(command_line, shell=True)
		return output_prefix + ".CNVcall.1000.vcf"

	def Manta_config(
		self, output_prefix):
		"""
		Generates Manta config

		Parameters
		----------

		output_prefix: str
			path and prefix for output files

		Returns
		-------

		str
			path to manta directory

		Raises
		------

		To do!

		"""

		command_line = "configManta.py --bam {} --reference {} --runDir {}_Manta".format(
		self.file_path, self.reference, output_prefix)
		subprocess.call(command_line, shell=True)
		return output_prefix + '_Manta/'


	def run_Manta(
		self, Manta_dir):

		"""
		Runs Manta on a bwa-mem aligned whole genome bam to detect SVs.

		Parameters
		----------

		Manta_dir: str
			path to Manta directory containing runWorkflow.py

		Returns
		-------

		str
			path to manta vcf

		Raises
		------

		To do!

		"""
		command_line = "python {} -m local -j 2".format(
				Manta_dir + "/runWorkflow.py")
		subprocess.call(command_line, shell=True)

		command_line = "gunzip {}".format(Manta_Dir + "/results/variants/diploid.vcf.gz")
		subprocess.call(command_line, shell=True)

		return Manta_dir + "/results/variants/diploid.vcf"

	def run_ERDS(self, ERDS_path, output_prefix, small_variant_vcf):
		"""
		Generates Manta config

		Parameters
		----------
		ERDS_path : str
			path to erds_pipeline.pl

		output_prefix: str
			path and prefix for output files

		Returns
		-------

		str
			path to ERDS folder

		Raises
		------

		To do!

		"""
		command_line = "{} -b {} -v {} -o {} -r {} ".format(
				ERDS_path, self.file_path, small_variant_vcf, output_prefix + ".ERDS", self.reference)
		subprocess.call(command_line, shell=True)
		return output_prefix + '.ERDS/'

	def make_variant_file(
		self,output_prefix, *args):
		"""
		Generates list of paths to vcf files generated by SV callers

		Parameters
		----------

		output_prefix: str
			path to output directory including sample name

		*args : str
			paths to SV vcf files

		Returns
		-------

		str
			path to list of variants

		Raises
		------

		To do!

		"""
		variant_list = open(output_prefix + "_variant_list.txt", 'w')

		vcfs = ""

		for arg in args:
			vcfs = vcfs + arg + '\n'

		variant_list.write(vcfs)
		variant_list.close()

		return output_prefix + "_variant_list.txt"


	def run_SURVIVOR(
		self,SURVIVOR_path, variant_list, max_dist, min_support, type, strand, estimate_distance, min_SV, output_prefix):
		"""
		Generates list of paths to vcf files generated by SV callers

		Parameters
		----------
		SURVIVOR_path: str
			path to SURVIVOR executable

		variant_list: str
			path to txt file of SV vcfs

		max_dist: int
			max distance between breakpoints

		min_support: int
			minimum number of supporting callers

		type: 0 or 1
			take the type into account (1=yes, 0=no)

		strand: 0 or 1
			take the strands of SVs into account (1=yes, 0=no)

		estimate_distance: 0 or 1
			estimate distance based on the size of SV (1=yes, 0=no)

		min_SV: int
			minimum size of SVs to be taken into account

		output_prefix: str
			path to output directory including sample name


		Returns
		-------

		str
			path to list of merged SURVIVOR variants

		Raises
		------

		To do!

		"""
		command_line = "{} merge {} {} {} {} {} {} {} {}".format(SURVIVOR_path, variant_list,
max_dist, min_support, type, strand, estimate_distance, min_support, output_prefix + ".SURVIVOR.SV.vcf")
		subprocess.call(command_line, shell=True)
		return output_prefix + ".SURVIVOR.SV.vcf"
