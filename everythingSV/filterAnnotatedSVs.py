import sys
import re
import argparse

def parse_args():
    """ 
    Parse command-line arguments

    Returns
    ----
    Parser argument namespace
    """
    parser = argparse.ArgumentParser(
            description="This script filters a set of annotated SVs output from everythingSV.py")
    parser.add_argument(
            "--annovar_output_with_overlap", type=str, required=True,
            help="Path to annovar output with overlap e.g.<SAMPLE>.DEL.withoverlap.header.annovar.hg19_multianno.txt")
    parser.add_argument(
            "--workdir", type=str, required=True,
            help="Path of working directory")
    parser.add_argument(
            "--output_prefix", type=str, required=True,
            help="Sample prefix")
    parser.add_argument(
            "--svtype", type=str, required=True,
            help="SV type")
    parser.add_argument(
            "--max_AC",type=int,  required=True,
            help="Variants with an allele count greater than max_AC will be filtered out if the reciprocal overlap with the sample SV is greater than or equal to 0.8")
    args = parser.parse_args()
    return args

def main():
	args = parse_args()
	SVs = open(args.annovar_output_with_overlap, 'r')
	filtered_SVs = open(args.workdir + "/" + args.output_prefix + "." + args.svtype +  ".filtered.withoverlap.header.annovar.hg19_multianno.txt", 'w')
	for line in SVs:
		cols = line.strip('\n').split('\t')
		CAUSESinhouse = cols[10]
		IMAGINEinhouse = cols[11]
		if cols[0] == "Chr":
			filtered_SVs.write("%s"%line)
		elif CAUSESinhouse == "":
			if IMAGINEinhouse == "":
				#no overlap with either database
				filtered_SVs.write("%s"%line)
			else:
				#overlap only with IMAGINE database
				IMAGINEinhouse =IMAGINEinhouse.split(',')
				AClist = []
				for record in IMAGINEinhouse: 
					record = record.split(':')
					AC = int(record[3])
					overlap = float(record[4])
					if overlap >= 0.8:
						AClist.append(AC)
					else:
						pass
				#IMAGINE count > max_AC
				if len(AClist) == 0:
					filtered_SVs.write("%s"%line)
				elif max(AClist) > args.max_AC:
					pass
				else:
					filtered_SVs.write("%s"%line)
		elif IMAGINEinhouse == "": 
			#overlap only with CAUSES database
			CAUSESinhouse =CAUSESinhouse.split(',')
			AClist = []
			for record in CAUSESinhouse: 
				record = record.split(':')
				AC = int(record[3])
				overlap = float(record[4])
				if overlap >= 0.8:
					AClist.append(AC)
				else:
					pass
			if len(AClist) == 0:
				filtered_SVs.write("%s"%line)
			elif max(AClist) > args.max_AC:
				pass
			else:
				filtered_SVs.write("%s"%line)
		else: 
			#overlap with both databases
			CAUSESinhouse =CAUSESinhouse.split(',')
			AClist = []
			for record in CAUSESinhouse: 
				record = record.split(':')
				AC = int(record[3])
				overlap = float(record[4])
				if overlap >= 0.8:
					AClist.append(AC)
				else:
					pass
			IMAGINEinhouse =IMAGINEinhouse.split(',')
			AClist = []
			for record in IMAGINEinhouse: 
				record = record.split(':')
				AC = int(record[3])
				overlap = float(record[4])
				if overlap >= 0.8:
					AClist.append(AC)
				else:
					pass
			if len(AClist) == 0:
					filtered_SVs.write("%s"%line)
			elif max(AClist) > args.max_AC:
				pass
			else:
				filtered_SVs.write("%s"%line)


main()
















