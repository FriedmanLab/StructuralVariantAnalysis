import sys
import re
import argparse
import bam
def parse_args():

	"""
	Parse command-line arguments

	Returns
	------

	Parser argument namespace
	"""

	parser = argparse.ArgumentParser(
			description="This script runs LUMPY, Manta, CNVnator, and ERDS"
			"on a bwa-mem aligned bam file." )
	parser.add_argument(
			"--bam_file", type=str, required=True,
			help="Path to bam file")
	parser.add_argument(
			"--reference", type=str, required=True,
			help="Path to reference fasta")
	parser.add_argument(
			"--chromosomes", type=str, required=True,
			help="Path to chromosome fastas")
	parser.add_argument(
		"--samtools_path", type=str, required=True,
		help="Path to Samtools executable")
	parser.add_argument(
		"--working_dir", type=str, required=True,
		help="Path to working directory")
	parser.add_argument(
		"--output_prefix", type=str, required=True,
		help="Output prefix for filenames, e.g. sample ID")
	#LUMPY-specific arguments
	parser.add_argument(
		"--lumpy_path", type=str, required=True,
		help="Path to LUMPY directory")
	parser.add_argument(
		"--lumpy", type=str, required=True,
		help="Path to LUMPY executable")
	parser.add_argument(
		"--read_length", type=str, required=True,
		help="Read length")
	parser.add_argument(
		"--min_mapq", type=str, required=True,
		help="Minimum mapping quality for running LUMPY")
	parser.add_argument(
		"--min_support", type=str, required=True,
		help="Minimum number of support reads for LUMPY call")
	#Manta-specific arguments
	parser.add_argument(
		"--manta_config", type=str, required=True,
		help="path to configManta.py (i.e. (MANTA_INSTALL_PATH)/bin/configManta.py)")
	#ERDS-specific arguments
	parser.add_argument(
		"--ERDS_path", type=str, required=True,
		help="Path to erds_pipeline.pl")
	parser.add_argument(
		"--small_variant_vcf", type=str, required=True,
		help="Path to GATK SNV/indel vcf")
	#SURVIVOR-specific arguments
	parser.add_argument(
		"--SURVIVOR_path", type=str, required=True,
		help="Path to SURVIVOR executable")
	parser.add_argument(
		"--SURVIVOR_max_dist", type=str, required=True,
		help="Max distance between breakpoints")
	parser.add_argument(
		"--SURVIVOR_min_support", type=str, required=True,
		help="Minimum number of supporting callers")
	parser.add_argument(
		"--SURVIVOR_type", type=str, required=True,
		help="take the type into account (1=yes, 0=no)")
	parser.add_argument(
		"--SURVIVOR_strand", type=str, required=True,
		help="take the strands of SVs into account (1=yes, 0=no)")
	parser.add_argument(
		"--SURVIVOR_estimate_distance", type=str, required=True,
		help="estimate distance based on the size of SV (1=yes, 0=no)")
	parser.add_argument(
		"--SURVIVOR_min_SV", type=str, required=True,
		help="minimum size of SVs to be taken into account")


	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	#instantiate bam object
	sample_bam = bam.bam_file(args.bam_file, args.reference, args.samtools_path)

	#Extract and sort discordant and split reads
	#sample_bam.extract_discordants(args.output_prefix)
	#sample_bam.extract_splitters(args.lumpy_path, args.output_prefix)
	discordants = args.working_dir + args.output_prefix + ".discordants.unsorted.bam"
	splitters = args.working_dir + args.output_prefix + ".splitters.unsorted.bam"
	splitters_sorted = args.output_prefix + ".splitters"
	discordants_sorted = args.output_prefix + ".discordants"
	#sample_bam.sort_bam(splitters, splitters_sorted)
	#sample_bam.sort_bam(discordants, discordants_sorted)
	insert = sample_bam.get_insert_size(args.lumpy_path, args.output_prefix)

	#Run LUMPY
	LUMPY = sample_bam.run_LUMPY(args.lumpy, args.output_prefix,
		sample_bam.get_insert_size(args.lumpy_path, args.output_prefix), args.read_length, args.min_mapq, args.min_support)

	#Run CNVnator
	CNVnator = sample_bam.run_CNVnator(args.chromosomes, args.output_prefix)

	#Generate Manta config file, then run Manta
	Manta_Dir = sample_bam.Manta_config(args.manta_config, args.output_prefix)
	Manta = sample_bam.run_Manta(Manta_Dir)

	#Run ERDS
	ERDS = sample_bam.run_ERDS(args.ERDS_path, args.output_prefix, args.small_variant_vcf)

	#Generate list of SV vcfs, then run SURVIVOR
	variant_list = sample_bam.make_variant_file(args.output_prefix, LUMPY, CNVnator, Manta, ERDS)

	SURVIVOR = sample_bam.run_SURVIVOR(args.SURVIVOR_path, variant_list,
	args.SURVIVOR_max_dist, args.SURVIVOR_min_support,
	args.SURVIVOR_type, args.SURVIVOR_strand, args.SURVIVOR_estimate_distance, args.SURVIVOR_min_SV, args.output_prefix)







main()
