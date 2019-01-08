# everythingSV

everythingSV comprises a pipeline for comprehensive structural variant (SV) analysis from whole genome sequences. The pipelines performs SV calling, merging, and annotation. This pipeline leverages split reads, discordant read pairs, and read depth information by utilizing four SV callers: [LUMPY](https://github.com/arq5x/lumpy-sv), [CNVnator](https://github.com/abyzovlab/CNVnator), [Manta](https://github.com/Illumina/manta), and [ERDS](http://people.duke.edu/~mz34/erds.htm). [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) is used to merge vcfs from the aforementioned callers. [Annovar](http://annovar.openbioinformatics.org/en/latest/) is used to annotate the SVs. The pipeline divides the SVs into 'small' and 'large' SVs by allowing the user the specify the minimum SV size for SURVIVOR merging (e.g. 50 bp for 'small' SVs and 50kb for 'large' SVs). The 'small' SVs include deletions, duplications, inversions, insertions, and unclassified breakends; this SV set represents variants not typically identified on chromosomal microarray. The 'large' SVs include deletions and duplications called by CNVnator and ERDS, and are meant to identify CNVs typically identified on chromosomal microarray. SVs are automatically annotated against RefSeq genes, and the user can provide additional custom bed files against which to annotate sample SVs, for instance databases of common CNVs or repetitive regions. 

Some of this pipeline was created in collaboration with Phillip Richmond.

## Software Requirements
* [python 2.7](https://www.python.org/)
  * [pandas](http://pandas.pydata.org/)
* [samtools](http://www.htslib.org/)
* [LUMPY](https://github.com/arq5x/lumpy-sv)
* [CNVnator](https://github.com/abyzovlab/CNVnator)
* [Manta](https://github.com/Illumina/manta)
* [ERDS](http://people.duke.edu/~mz34/erds.htm)
* [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)

## Installation

'''
git clone https://github.com/FriedmanLab/StructuralVariantAnalysis

'''
