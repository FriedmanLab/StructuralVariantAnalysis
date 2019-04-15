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
            description="This script filters a set of annotated SVs output from everythingSV.py based on overlap with inhouse databases")
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
            help="Variants that have a maximum allele count greater than max_AC, or overlap with more than max_AC number of different variants, are filtered out. e.g. say your max_AC is 5 (1) if AC are 6,5,4,3,2,1 then this would be filtered because AC1 = 6, (2) 5,4,3,2,1,1 would also be filtered because there are 6 'different' AC ")
    args = parser.parse_args()
    return args

def main():
	args = parse_args()
	SVs = open(args.annovar_output_with_overlap, 'r')
	filtered_SVs = open(args.workdir + "/" + args.output_prefix + "." + args.svtype +  ".filtered.withoverlap.header.annovar.hg19_multianno.txt", 'w')
	for line in SVs:
		cols = line.strip('\n').split('\t')
		CAUSESinhouse = cols[10].split(',')
		IMAGINEinhouse = cols[11].split(',')
		if cols[0] == "Chr":
			filtered_SVs.write("%s"%line)
		elif CAUSESinhouse[0] == "":
			if IMAGINEinhouse[0] == "":
				#no overlap with either database
				filtered_SVs.write("%s"%line)
			else:
				#overlap only with IMAGINE database
				IMAGINEinhouse = map(int, IMAGINEinhouse)
				if max(IMAGINEinhouse) > args.max_AC:
					pass
				elif len(IMAGINEinhouse) > args.max_AC:
					pass
				else:
					filtered_SVs.write("%s"%line)
		elif IMAGINEinhouse[0] == "": 
			#overlap only with CAUSES database
			CAUSESinhouse = map(int, CAUSESinhouse)
			if max(CAUSESinhouse) > int(args.max_AC):
				pass
			elif len(CAUSESinhouse) > int(args.max_AC):
				pass
			else:
				filtered_SVs.write("%s"%line)
		else: 
			CAUSESinhouse = map(int, CAUSESinhouse)
			IMAGINEinhouse = map(int, IMAGINEinhouse)
			#overlap with both databases
			if max(CAUSESinhouse) > int(args.max_AC):
				pass
			elif len(CAUSESinhouse) > int(args.max_AC):
				pass
			else:
				if max(IMAGINEinhouse) > int(args.max_AC):

					pass
				elif len(IMAGINEinhouse) > int(args.max_AC):
					pass
				else:
					filtered_SVs.write("%s"%line)


main()
















