import sys
import argparse

#convert refGene file into bed format 


def parse_args():
	"""
	Parse command-line arguments

	Returns
	------
	Parser argument namespace
	"""

	parser = argparse.ArgumentParser(description = "This script takes a list of variants output from VarSeq and parses into a bed file with desired annotations")

	parser.add_argument("--variants", type = str, required = True, help = "path to VarSeq variant file in tsv format")

	parser.add_argument("--coordinates", type = int, required = True, help = "column number for variant coordinates")

	parser.add_argument("--annotations", required = True, nargs='+', help = "A list where the column number of the desired annotation is the first element,and the second element is the name of the annotation e.g.'6,CADD_score'.")
	parser.add_argument("--output", type = str, required = True, help = "Name of output file")

	args = parser.parse_args()
	return args


def main():
	args = parse_args()
	#variant coordinates
	coordinates = int(args.coordinates) - 1

	variants = open(args.variants, 'r')
	variants_out = open(args.output, 'w')


	for variant in variants:
		var_cols = variant.strip('\n').split('\t')
		if var_cols[0] == 'Gene Name':
			pass
		else:
			#chromsome and position
			var_coord = var_cols[coordinates]
			chr = var_coord.split(':')[0]
			start = var_coord.split(':')[1]
			end = str(int(start) + 1)
			#grab specified annotations from variant row
			anno = []
			for col in args.annotations:
				#column index in variant file
				col_num = int(col.split(',')[0])-1
				#column title in variant file
				col_title = col.split(',')[1]
				col_value = var_cols[col_num]
				col_value = col_title + ":" + col_value
				anno.append(col_value)
			new_cols = [chr, start, end]
			var_coord = "coordinates:" + var_coord
			anno.append(var_coord)
			anno = ";".join(anno)
			new_cols.append(anno)
			newline = "\t".join(new_cols)
			variants_out.write("%s\n"%newline)








main()