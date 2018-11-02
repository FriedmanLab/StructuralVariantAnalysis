
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
        DEL = workdir + output_prefix + ".DEL"
        DEL = open(DEL, 'w')
        DUP = workdir + output_prefix + ".DUP"
        DUP = open(DUP, 'w')
        INV = workdir + output_prefix + ".INV"
        INV = open(INV, 'w')
        INS = workdir + output_prefix + ".INS"
        INS = open(INS, 'w')
        BND = workdir + output_prefix + ".BND"
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
            A list where the file name of regions of interest in bed format is the first element, an integer specifiying the column of interest for annotation is the second element, and the minimum fraction of overlap of sample SV with region of interest is the third element (e.g. "gene_promoters.bed, 4, 0.9")

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
            output_prefix + '.annovar.txt', protocols, filenames, operations, bed_args)

        print command_line

        subprocess.call(command_line, shell=True)

        return output_prefix + ".annovar.txt.hg19_multianno.txt"

