import sys
import re
import argparse
import subprocess
import bam
import SURVIVORvcf
import annotateSURVIVOR

def parse_args():

	"""
	Parse command-line arguments

	Returns
	------

	Parser argument namespace
	"""

	parser = argparse.ArgumentParser(
		description="This script runs LUMPY, Manta, CNVnator, and ERDS"
		"on a bwa-mem aligned bam file. It merges vcfs from these tools using SURVIVOR, and annotates SVs using Annovar." )
	parser.add_argument(
		"--bam_file", type=str, required=True,
		help="Path to bam file. Required for SV calling and annotation pipelines.")
	parser.add_argument(
		"--reference", type=str, required=True,
		help="Path to reference fasta. Required for SV calling and annotation pipelines.")
	parser.add_argument(
		"--samtools_path", type=str, required=True,
		help="Path to Samtools executable. Required for SV calling and annotation pipelines.")
	parser.add_argument(
		"--working_dir", type=str, required=True,
		help="Path to working directory. Required for SV calling and annotation pipelines. ")
	parser.add_argument(
		"--output_prefix", type=str, required=True,
		help="Output prefix for filenames, e.g. sample ID. Required for SV calling and annotation pipelines.")
	
	#LUMPY-specific arguments
	parser.add_argument(
		"--lumpy_path", type=str, required=False,
		help="Path to LUMPY directory. Required for SV calling pipeline.")
	parser.add_argument(
		"--lumpy", type=str, required=False,
		help="Path to LUMPY executable. Required for SV calling pipeline.")
	parser.add_argument(
		"--read_length", type=str, required=False,
		help="Read length. Required for SV calling pipeline.")
	parser.add_argument(
		"--min_mapq", type=str, required=False,
		help="Minimum mapping quality for running LUMPY. Required for SV calling pipeline.")
	parser.add_argument(
		"--min_support", type=str, required=False,
		help="Minimum number of support reads for LUMPY call. Required for SV calling pipeline.")

	#CNVnator-specific arguments
	parser.add_argument(
		"--chromosomes", type=str, required=False,
		help="Path to chromosome fastas. Required for SV calling pipeline.")
	parser.add_argument(
		"--CNVnator", type=str, required=False,
		help="Path to CNVnator exectuable. Required for SV calling pipeline.")
	parser.add_argument(
		"--cnvnator2VCF", type=str, required=False,
		help="Path to cnvnator2VCF.pl script. Required for SV calling pipeline.")
	parser.add_argument(
		"--bcftools", type=str, required=False,
		help="Path to bcftools exectuable. Required for SV calling pipeline.")

	#Manta-specific arguments
	parser.add_argument(
		"--manta_config", type=str, required=False,
		help="path to configManta.py (i.e. (MANTA_INSTALL_PATH)/bin/configManta.py). Required for SV calling pipeline.")

	#ERDS-specific arguments
	parser.add_argument(
		"--ERDS_path", type=str, required=False,
		help="Path to erds_pipeline.pl. Required for SV calling pipeline.")
	parser.add_argument(
		"--small_variant_vcf", type=str, required=False,
		help="Path to GATK SNV/indel vcf. Required for SV calling pipeline.")

	#SURVIVOR-specific arguments
	parser.add_argument(
		"--SURVIVOR_path", type=str, required=False,
		help="Path to SURVIVOR executable. Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_max_dist_small", type=str, required=False,
		help="Max distance between breakpoints for 'small' SVs. Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_max_dist_large", type=str, required=False,
		help="Max distance between breakpoints for 'large' SVs. Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_min_support", type=str, required=False,
		help="Minimum number of supporting callers. Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_type", type=str, required=False,
		help="take the type into account (1=yes, 0=no). Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_strand", type=str, required=False,
		help="take the strands of SVs into account (1=yes, 0=no). Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_estimate_distance", type=str, required=False,
		help="estimate distance based on the size of SV (1=yes, 0=no). Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_min_SV_small", type=str, required=False,
		help="minimum size of 'small' SVs to be taken into account. Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_min_SV_large", type=str, required=False,
		help="minimum size of 'large' SVs to be taken into account. Required for annotation pipeline.")
	parser.add_argument(
		"--SURVIVOR_additional_vcfs", nargs='+', required=False,
		help="additional SV vcfs to input to SURVIVOR. Optional for annotation pipeline.")

	#Options for enabling or disabling parts of pipeline 
	parser.add_argument(
		"--SV_calling_only", action="store_true", default=False,
		help="If true, perform SV valling with annotation")
	parser.add_argument(
		"--annotation_only", action="store_true", default=False,
		help="If true, perform annotation without SV calling")
	parser.add_argument(
		"--LUMPY_vcf", type=str, required=False,
		help="Path to LUMPY vcf; provide this if --annotation_only True")
	parser.add_argument(
		"--CNVnator_vcf", type=str, required=False,
		help="Path to CNVnator vcf; provide this if --annotation_only True")
	parser.add_argument(
		"--Manta_vcf", type=str, required=False,
		help="Path to Manta vcf; provide this if --annotation_only True")
	parser.add_argument(
		"--ERDS_vcf", type=str, required=False,
		help="Path to ERDS vcf; provide this if --annotation_only True")

	#Annovar-specific arguments
	parser.add_argument(
		"--table_annovar", type=str, required=False,
		help="path to Annovar table_annovar.pl script. Required for annotation pipeline.")
	parser.add_argument(
		"--humandb", type=str, required=False,
		help="path to Annovar humandb/ folder containing annotations of interest. Required for annotation pipeline.")
	parser.add_argument(
		"--annovar_reference", type=str, required=False,
		help="hg19 or hg38. Required for annotation pipeline.")
	parser.add_argument(
		"--deletion_bedfiles", required=False, nargs='+', help="A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, the minimum fraction of overlap of the sample SV with region of interest is the third element, the fourth element is TRUE if it is desired to have the percentage overlap of the database annotation calculated, otherwise FALSE, and the fifth element is the name of the annotation  (e.g. 'gene_promoters.bed, 4, 0.9, TRUE, Gene_promoters'. Required for annotation pipeline.")
	parser.add_argument(
		"--duplication_bedfiles", required=False, nargs='+', help="A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, the minimum fraction of overlap of the sample SV with region of interest is the third element,the fourth element is TRUE if it is desired to have the percentage overlap of the database annotation calculated, otherwise FALSE, and the fifth element is the name of the annotation  (e.g. 'gene_promoters.bed, 4, 0.9, TRUE, Gene_promoters'. Required for annotation pipeline.")
	parser.add_argument(
		"--inversion_bedfiles", required=False, nargs='+', help="A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, the minimum fraction of overlap of the sample SV with region of interest is the third element,the fourth element is TRUE if it is desired to have the percentage overlap of the database annotation calculated, otherwise FALSE, and the fifth element is the name of the annotation  (e.g. 'gene_promoters.bed, 4, 0.9, TRUE, Gene_promoters'. Required for annotation pipeline.")
	parser.add_argument(
		"--insertion_bedfiles", required=False, nargs='+', help="A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, the minimum fraction of overlap of the sample SV with region of interest is the third element,the fourth element is TRUE if it is desired to have the percentage overlap of the database annotation calculated, otherwise FALSE, and the fifth element is the name of the annotation  (e.g. 'gene_promoters.bed, 4, 0.9, TRUE, Gene_promoters'. Required for annotation pipeline.")
	parser.add_argument(
		"--breakend_bedfiles", required=False, nargs='+', help="A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, the minimum fraction of overlap of the sample SV with region of interest is the third element, the fourth element is TRUE if it is desired to have the percentage overlap of the database annotation calculated, otherwise FALSE, and the fifth element is the name of the annotation  (e.g. 'gene_promoters.bed, 4, 0.9, TRUE, Gene_promoters'. Required for annotation pipeline.")



	args = parser.parse_args()
	return args

def run_SV_pipeline(
	bam_file, reference, samtools_path, working_dir, output_prefix, lumpy_path, lumpy, read_length, min_mapq, min_support, CNVnator, cnvnator2VCF, bcftools, chromosomes, manta_config, ERDS_path, small_variant_vcf):
	"""
	Runs SV callers: LUMPY, CNVnator, Manta, ERDS.

	Parameters
	----------

	bam_file: str
		Path to bamfile
	reference: str
		Path to reference genome fasta 
	samtools_path: str 
		Path to samtools executable
	working_dir: str
		Path to working directory
	output_prefix: str
		Output prefix for filenames, e.g. sample ID
	lumpy_path: str 
		Path to LUMPY directory
	lumpy: str
		Path to LUMPY executable
	read_length: int
		Read length 
	min_mapq: int
		Minimum mapping quality for running LUMPY
	min_support: int
		Minimum number of support reads for LUMPY call
	CNVnator: str
		Path to CNVnator exectuable
	cnvnator2VCF: str
		Path to cnvnator2VCF.pl script
	bcftools: str
		Path to bcftools exectuable
	chromosomes: str 
		Path to chromosome fastas
	manta_config: str
		Path to configManta.py (i.e. (MANTA_INSTALL_PATH)/bin/configManta.py
	ERDS_path: str
		Path to ERDS exectuable
	small_variant_vcf: str 
		Path to GATK SNV/indel vcf

	Returns
	-------
	list
		list of SV vcfs where the first element is the LUMPY vcf, second is the CNVnator vcf, third is the Manta vcf, and fourth is the ERDS vcf
	"""

	#instantiate bam object
	sample_bam = bam.bam_file(bam_file, reference, samtools_path)

	#Extract and sort discordant and split reads
	sample_bam.extract_discordants(working_dir, output_prefix)
	sample_bam.extract_splitters(lumpy_path, working_dir, output_prefix)
	discordants = working_dir + '/' + output_prefix + ".discordants.unsorted.bam"
	splitters = working_dir + '/' + output_prefix + ".splitters.unsorted.bam"
	splitters_sorted = output_prefix + ".splitters"
	discordants_sorted = output_prefix + ".discordants"
	sample_bam.sort_bam(splitters, working_dir, splitters_sorted)
	sample_bam.sort_bam(discordants, working_dir, discordants_sorted)

	#Run LUMPY
	LUMPY = sample_bam.run_LUMPY(lumpy, working_dir, output_prefix,
		sample_bam.get_insert_size(lumpy_path, working_dir, output_prefix), read_length, min_mapq, min_support)

	#Run CNVnator
	CNVnator = sample_bam.run_CNVnator(CNVnator, cnvnator2VCF, bcftools, chromosomes, working_dir, output_prefix)

	#Generate Manta config file, then run Manta
	Manta_Dir = sample_bam.Manta_config(manta_config, working_dir, output_prefix)
	Manta = sample_bam.run_Manta(Manta_Dir)

	#Run ERDS
	ERDS = sample_bam.run_ERDS(ERDS_path, working_dir, output_prefix, small_variant_vcf)

	SV_list = [LUMPY, CNVnator, Manta, ERDS]
	
	return SV_list

def annotate_variants(
	working_dir, output_prefix, bam_file, reference, samtools_path, LUMPY, CNVnator, Manta, ERDS, SURVIVOR_additional_vcfs, SURVIVOR_path, SURVIVOR_max_dist_small,SURVIVOR_max_dist_large, SURVIVOR_min_support, SURVIVOR_type, SURVIVOR_strand, SURVIVOR_estimate_distance, SURVIVOR_min_SV_small, SURVIVOR_min_SV_large, table_annovar, humandb, annovar_reference, deletion_bedfiles, duplication_bedfiles, inversion_bedfiles, insertion_bedfiles, breakend_bedfiles):

	"""
	Merges SV vcfs via SURVIVOR, and annotates variants

	Paramters
	---------
	working_dir: str
		Path to working directory
	output_prefix: str
		Output prefix for filenames, e.g. sample ID
	bam_file: str
		Path to bam file
	reference: str
		Path to reference fasta
	samtools_path: str
		Path to samtools executable
	LUMPY: str
		Path to LUMPY vcf
	CNVnator: str 
		Path to CNVnator vcf
	Manta: str
		Path to Manta vcf
	ERDS: str
		Path to ERDS vcf
	SURVIVOR_additional_vcfs: str 
		Additional SV vcfs to input to SURVIVOR
	SURVIVOR_path: str
		Path to SURVIVOR exectuable
	SURVIVOR_max_dist_small: int
		Max distance between breakpoints for 'small' SVs
	SURVIVOR_max_dist_large: int
		Max distance between breakpoints for 'large' SVs
	SURVIVOR_min_support: int
		Minimum number of supporting callers
	SURVIVOR_type: int
		Take the type into account (1=yes, 0=no)
	SURVIVOR_strand: int
		Take the strands of SVs into account (1=yes, 0=no)
	SURVIVOR_estimate_distance: int
		Estimate distance based on the size of SV (1=yes, 0=no)
	SURVIVOR_min_SV_small: int
		Minimum size of 'small' SVs to be taken into account
	SURVIVOR_min_SV_large: int
		Minimum size of 'large' SVs to be taken into account
	table_annovar: str
		Path to Annovar table_annovar.pl script
	humandb: str
		Path to Annovar humandb/ folder containing annotations of interest
	annovar_reference: str 
		hg19 or hg38
	deletion_bedfiles: list
		A list of annotations against which sample SVs are annotated
	duplication_bedfiles: list
		A list of annotations against which sample SVs are annotated
	inversion_bedfiles: list
		A list of annotations against which sample SVs are annotated
	insertion_bedfiles: list
		A list of annotations against which sample SVs are annotated
	breakend_bedfiles: list
		A list of annotations against which sample SVs are annotated


	"""

	#instantiate bam object
	sample_bam = bam.bam_file(bam_file, reference, samtools_path)

	#ANNOTATE SMALL SVS
	#Generate list of SV vcfs, then run SURVIVOR
	if  SURVIVOR_additional_vcfs != None:
		variant_list = sample_bam.make_variant_file(working_dir, output_prefix, LUMPY, Manta, SURVIVOR_additional_vcfs)
	else:
		variant_list = sample_bam.make_variant_file(working_dir, output_prefix, LUMPY, Manta)
	
	SURVIVOR = sample_bam.run_SURVIVOR(SURVIVOR_path, variant_list,
	SURVIVOR_max_dist_small, SURVIVOR_min_support,
	SURVIVOR_type, SURVIVOR_strand, SURVIVOR_estimate_distance, SURVIVOR_min_SV_small, working_dir, output_prefix, 'small')

	#Convert SURIVIVOR vcf to avinput
	SURVIVOR = SURVIVORvcf.SURVIVOR_vcf(SURVIVOR)
	avinput = SURVIVOR.SURVIVOR_to_avinput(output_prefix, working_dir, 'small')
	
	#Separate SVs by type (ie DEL, DUP, INS, INV, BND)
	annotateSURVIVOR.split_by_SV_type(avinput, output_prefix, working_dir, 'small')
	DEL = working_dir + "/" +  output_prefix + ".DEL." + 'small'
	DUP = working_dir +  "/" + output_prefix + ".DUP." + 'small'
	INV = working_dir +  "/" + output_prefix + ".INV." + 'small'
	INS = working_dir +  "/" + output_prefix + ".INS." + 'small'
	BND = working_dir +  "/" + output_prefix + ".BND." + 'small'

	#Annotate variants
	annotated_DEL = annotateSURVIVOR.annotate_avinput(DEL, table_annovar, humandb, annovar_reference, output_prefix + ".DEL", working_dir, 'small', deletion_bedfiles) 

	annotated_DUP = annotateSURVIVOR.annotate_avinput(DUP, table_annovar, humandb, annovar_reference, output_prefix + ".DUP", working_dir,'small', duplication_bedfiles) 

	annotated_INV = annotateSURVIVOR.annotate_avinput(INV, table_annovar, humandb, annovar_reference, output_prefix + ".INV", working_dir, 'small', inversion_bedfiles) 

	annotated_INS = annotateSURVIVOR.annotate_avinput(INS, table_annovar, humandb, annovar_reference, output_prefix + ".INS", working_dir, 'small',insertion_bedfiles) 

	annotated_BND = annotateSURVIVOR.annotate_avinput(BND, table_annovar, humandb, annovar_reference, output_prefix + ".BND", working_dir, 'small', breakend_bedfiles) 

	#Calculate reciprocal overlap for specified bed files
	annotated_DEL_overlap = annotateSURVIVOR.calculate_overlap(output_prefix, working_dir, annotated_DEL, "DEL", 'small',  deletion_bedfiles)

	annotated_DUP_overlap = annotateSURVIVOR.calculate_overlap(output_prefix, working_dir,annotated_DUP, "DUP", 'small', duplication_bedfiles)

	annotated_INV_overlap = annotateSURVIVOR.calculate_overlap(output_prefix, working_dir,annotated_INV, "INV", 'small',inversion_bedfiles)

	annotated_INS_overlap = annotateSURVIVOR.calculate_overlap(output_prefix, working_dir,annotated_INS, "INS", 'small',  insertion_bedfiles)

	annotated_BND_overlap = annotateSURVIVOR.calculate_overlap(output_prefix, working_dir,annotated_BND, "BND",'small',  breakend_bedfiles)

	#Add headers
	#Deletion
	DEL_header = annotateSURVIVOR.add_header(output_prefix, working_dir, "DEL.small",  deletion_bedfiles)

	command_line = "cat {} {} > {} ".format(DEL_header, annotated_DEL_overlap, working_dir + "/" +  output_prefix + ".DEL.small" + ".withoverlap.header.annovar.hg19_multianno.txt")
	subprocess.call(command_line, shell=True) 

	#Duplication
	DUP_header = annotateSURVIVOR.add_header(output_prefix, working_dir, "DUP.small",  duplication_bedfiles)

	command_line = "cat {} {} > {} ".format(DUP_header, annotated_DUP_overlap, working_dir + "/" + output_prefix + ".DUP.small" + ".withoverlap.header.annovar.hg19_multianno.txt")
	subprocess.call(command_line, shell=True)

	#Inversion
	INV_header = annotateSURVIVOR.add_header(output_prefix, working_dir, "INV.small",  inversion_bedfiles)

	command_line = "cat {} {} > {} ".format(INV_header, annotated_INV_overlap, working_dir + "/" + output_prefix + ".INV.small" + ".withoverlap.header.annovar.hg19_multianno.txt")
	subprocess.call(command_line, shell=True)

	#Insertion
	INS_header = annotateSURVIVOR.add_header(output_prefix, working_dir, "INS.small",  insertion_bedfiles)

	command_line = "cat {} {} > {} ".format(INS_header, annotated_INS_overlap, working_dir + "/" + output_prefix + ".INS.small" + ".withoverlap.header.annovar.hg19_multianno.txt")
	subprocess.call(command_line, shell=True)

	#Breakend
	BND_header = annotateSURVIVOR.add_header(output_prefix, working_dir, "BND.small",  breakend_bedfiles)

	command_line = "cat {} {} > {} ".format(BND_header, annotated_BND_overlap, working_dir + "/" + output_prefix + ".BND.small" + ".withoverlap.header.annovar.hg19_multianno.txt")
	subprocess.call(command_line, shell=True)

	#ANNOTATE LARGE CNVS
	large_variant_list = sample_bam.make_variant_file(working_dir, output_prefix, CNVnator, ERDS)
	
	SURVIVORlarge = sample_bam.run_SURVIVOR(SURVIVOR_path, large_variant_list,
	SURVIVOR_max_dist_large, SURVIVOR_min_support,
	SURVIVOR_type, SURVIVOR_strand, SURVIVOR_estimate_distance, SURVIVOR_min_SV_large, working_dir, output_prefix, 'large')
	
	#Convert SURIVIVOR vcf to avinput
	SURVIVORlarge = SURVIVORvcf.SURVIVOR_vcf(SURVIVORlarge)
	avinput_large = SURVIVORlarge.SURVIVOR_to_avinput(output_prefix, working_dir, 'large')

	#Separate SVs by type (DEL, DUP)
	annotateSURVIVOR.split_by_SV_type(avinput_large, output_prefix, working_dir, 'large')
	DEL_large = working_dir + "/" +  output_prefix + ".DEL." + 'large'
	DUP_large = working_dir +  "/" + output_prefix + ".DUP." + 'large'

	#Annotate variants
	annotated_DEL_large = annotateSURVIVOR.annotate_avinput(DEL_large, table_annovar, humandb, annovar_reference, output_prefix + ".DEL", working_dir, 'large', deletion_bedfiles) 

	annotated_DUP_large = annotateSURVIVOR.annotate_avinput(DUP_large, table_annovar, humandb, annovar_reference, output_prefix + ".DUP", working_dir, 'large', duplication_bedfiles) 

	#Calculate reciprocal overlap for specified bed files
	annotated_DEL_overlap_large = annotateSURVIVOR.calculate_overlap(output_prefix, working_dir, annotated_DEL_large, "DEL", 'large',  deletion_bedfiles)

	annotated_DUP_overlap_large = annotateSURVIVOR.calculate_overlap(output_prefix, working_dir, annotated_DUP_large, "DUP", 'large',  duplication_bedfiles)

	#Add headers
	#Deletion
	DEL_header_large = annotateSURVIVOR.add_header(output_prefix, working_dir, "DEL.large",  deletion_bedfiles)

	command_line = "cat {} {} > {} ".format(DEL_header_large, annotated_DEL_overlap_large, working_dir + "/" +  output_prefix + ".DEL.large" + ".withoverlap.header.annovar.hg19_multianno.txt")
	subprocess.call(command_line, shell=True) 

	DUP_header_large = annotateSURVIVOR.add_header(output_prefix, working_dir, "DUP.large",  duplication_bedfiles)

	command_line = "cat {} {} > {} ".format(DUP_header_large, annotated_DUP_overlap_large, working_dir + "/" +  output_prefix + ".DUP.large" + ".withoverlap.header.annovar.hg19_multianno.txt")
	subprocess.call(command_line, shell=True) 


def main():
	args = parse_args()

	#SV calling
	if args.SV_calling_only == False: 
		if args.annotation_only == False:
			print "Running SV and Annotation Pipelines..."
			#Run SV calling and annotation pipelines
			if all(v is not None for v in [args.bam_file, args.reference, args.samtools_path, args.working_dir, args.output_prefix, args.lumpy_path, args.lumpy, args.read_length, args.min_mapq, args.min_support, args.CNVnator,args.cnvnator2VCF, args.bcftools, args.chromosomes, args.manta_config, args.ERDS_path, args.small_variant_vcf]):
				#create working directory where results will be written
				command_line = "mkdir {}".format(args.working_dir)
				subprocess.call(command_line, shell=True)
				
				SVs = run_SV_pipeline(
					bam_file = args.bam_file,
					reference = args.reference, 
					samtools_path = args.samtools_path,
					working_dir = args.working_dir,
					output_prefix = args.output_prefix,
					lumpy_path = args.lumpy_path,
					lumpy = args.lumpy,
					read_length = args.read_length,
					min_mapq = args.min_mapq,
					min_support = args.min_support,
					CNVnator = args.CNVnator,
					cnvnator2VCF = args.cnvnator2VCF,
					bcftools = args.bcftools,
					chromosomes = args.chromosomes,
					manta_config = args.manta_config,
					ERDS_path = args.ERDS_path,
					small_variant_vcf = args.small_variant_vcf)

				annotate_variants(
					working_dir = args.working_dir,
					output_prefix = args.output_prefix,
					bam_file = args.bam_file,
					reference = args.reference,
					samtools_path = args.samtools_path,
					LUMPY = SVs[0],
					CNVnator = SVs[1],
					Manta = SVs[2],
					ERDS = SVs[3],
					SURVIVOR_additional_vcfs = args.SURVIVOR_additional_vcfs,
					SURVIVOR_path = args.SURVIVOR_path,
					SURVIVOR_max_dist_small = args.SURVIVOR_max_dist_small,
					SURVIVOR_max_dist_large = args.SURVIVOR_max_dist_large, 
					SURVIVOR_min_support = args.SURVIVOR_min_support,
					SURVIVOR_type = args.SURVIVOR_type,
					SURVIVOR_strand = args.SURVIVOR_strand,
					SURVIVOR_estimate_distance = args.
					SURVIVOR_estimate_distance,
					SURVIVOR_min_SV_small = args.SURVIVOR_min_SV_small, 
					SURVIVOR_min_SV_large = args.SURVIVOR_min_SV_large, 
					table_annovar = args.table_annovar,
					humandb = args.humandb,
					annovar_reference = args.annovar_reference,
					deletion_bedfiles = args.deletion_bedfiles,
					duplication_bedfiles = args.duplication_bedfiles,
					inversion_bedfiles = args.inversion_bedfiles, 
					insertion_bedfiles = args.insertion_bedfiles,
					breakend_bedfiles = args.breakend_bedfiles)
			else:
				sys.exit("Missing one or more arguments. See 'python everythingSV.py -h' for details on arguments.")
		else:
			#Run annotation pipeline only
			print "Running Annotation Pipeline..."
			if all(v is not None for v in [args.bam_file, args.reference, args.samtools_path, args.working_dir, args.output_prefix, args.LUMPY_vcf, args.CNVnator_vcf, args.Manta_vcf, args.ERDS_vcf, args.SURVIVOR_path, args.SURVIVOR_max_dist_small, args.SURVIVOR_max_dist_large, args.SURVIVOR_min_support, args.SURVIVOR_type, args.SURVIVOR_strand, args.SURVIVOR_estimate_distance, args.SURVIVOR_min_SV_small, args.SURVIVOR_min_SV_large, args.table_annovar, args.humandb, args.annovar_reference, args.deletion_bedfiles, args.duplication_bedfiles, args.inversion_bedfiles,args.insertion_bedfiles, args.breakend_bedfiles]):
					annotate_variants(
						working_dir = args.working_dir,
						output_prefix = args.output_prefix,
						bam_file = args.bam_file,
						reference = args.reference,
						samtools_path = args.samtools_path,
						LUMPY = args.LUMPY_vcf,
						CNVnator = args.CNVnator_vcf,
						Manta = args.Manta_vcf,
						ERDS = args.ERDS_vcf,
						SURVIVOR_additional_vcfs = args.
						SURVIVOR_additional_vcfs,
						SURVIVOR_path = args.SURVIVOR_path,
						SURVIVOR_max_dist_small = args.SURVIVOR_max_dist_small,
						SURVIVOR_max_dist_large = args.SURVIVOR_max_dist_large,
						 SURVIVOR_min_support = args.SURVIVOR_min_support,
						SURVIVOR_type = args.SURVIVOR_type,
						SURVIVOR_strand = args.SURVIVOR_strand,
						SURVIVOR_estimate_distance = args.
						SURVIVOR_estimate_distance,
						SURVIVOR_min_SV_small = args.SURVIVOR_min_SV_small, 
						SURVIVOR_min_SV_large = args.SURVIVOR_min_SV_large, 
						table_annovar = args.table_annovar,
						humandb = args.humandb,
						annovar_reference = args.annovar_reference,
						deletion_bedfiles = args.deletion_bedfiles,
						duplication_bedfiles = args.duplication_bedfiles,
						inversion_bedfiles = args.inversion_bedfiles, 
						insertion_bedfiles = args.insertion_bedfiles,
						breakend_bedfiles = args.breakend_bedfiles)
			else:
				sys.exit("Missing one or more arguments. See 'python everythingSV.py -h' for details on arguments.")

	else:
		#Run SV pipeline only
		print "Running SV Pipeline..."
		if all(v is not None for v in [args.bam_file, args.reference, args.samtools_path, args.working_dir, args.output_prefix, args.lumpy_path, args.lumpy, args.read_length, args.min_mapq, args.min_support, args.CNVnator, args.cnvnator2VCF, args.bcftools, args.chromosomes, args.manta_config, args.ERDS_path, args.small_variant_vcf]):
				#create working directory where results will be written
				command_line = "mkdir {}".format(args.working_dir)
				subprocess.call(command_line, shell=True)
				SVs = run_SV_pipeline(
					bam_file = args.bam_file,
					reference = args.reference, 
					samtools_path = args.samtools_path,
					working_dir = args.working_dir,
					output_prefix = args.output_prefix,
					lumpy_path = args.lumpy_path,
					lumpy = args.lumpy,
					read_length = args.read_length,
					min_mapq = args.min_mapq,
					min_support = args.min_support,
					CNVnator = args.CNVnator,
					cnvnator2VCF = args.cnvnator2VCF,
					bcftools = args.bcftools,
					chromosomes = args.chromosomes,
					manta_config = args.manta_config,
					ERDS_path = args.ERDS_path,
					small_variant_vcf = args.small_variant_vcf)
		else:
			sys.exit("Missing one or more arguments. See 'python everythingSV.py -h' for details on arguments.")



main()
