# everythingSV

everythingSV comprises a pipeline for comprehensive structural variant (SV) analysis from whole genome sequences. The pipelines performs SV calling, merging, and annotation. This pipeline leverages split reads, discordant read pairs, and read depth information by utilizing four SV callers: [LUMPY](https://github.com/arq5x/lumpy-sv), [CNVnator](https://github.com/abyzovlab/CNVnator), [Manta](https://github.com/Illumina/manta), and [ERDS](http://people.duke.edu/~mz34/erds.htm). [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) is used to merge vcfs from the aforementioned callers. [Annovar](http://annovar.openbioinformatics.org/en/latest/) is used to annotate the SVs. The pipeline divides the SVs into 'small' and 'large' SVs by allowing the user the specify the minimum SV size for SURVIVOR merging (e.g. 50 bp for 'small' SVs and 50kb for 'large' SVs). The 'small' SVs include deletions, duplications, inversions, insertions, and unclassified breakends; this SV set represents variants not typically identified on chromosomal microarray. The 'large' SVs include deletions and duplications called by CNVnator and ERDS, and are meant to identify CNVs typically identified on chromosomal microarray. SVs are automatically annotated against RefSeq genes, and the user can provide additional custom bed files against which to annotate sample SVs, for instance databases of common CNVs or repetitive regions. 

Some of this pipeline was created in collaboration with Phillip Richmond.

## Software Requirements

everythingSV requires the tools listed below already be installed:

* [python 2.7](https://www.python.org/)
  * [pandas](http://pandas.pydata.org/)
* [samtools](http://www.htslib.org/)
* [LUMPY](https://github.com/arq5x/lumpy-sv)
* [CNVnator](https://github.com/abyzovlab/CNVnator)
* [Manta](https://github.com/Illumina/manta)
* [ERDS](http://people.duke.edu/~mz34/erds.htm)
* [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
* [Annovar](http://annovar.openbioinformatics.org/en/latest/)

The paths to the LUMPY, CNVnator, Manta, ERDS, SURVIVOR, and annovar executables will be used as arguments in everythingSV. 



## Installation
```
git clone https://github.com/FriedmanLab/StructuralVariantAnalysis
```

## everythingSV Pipeline

The pipeline consists of two modules, the SV calling module and the SV annotation module. If the ```--sv_calling_only``` flag is used, SVs will be called annotation. If the ```--annotation_only``` flag is used, SVs are merged with SURVIVOR and annotated. If neither flag is specified, SVs are called, merged, and annotated. 

### SV calling only

```

python everythingSV.py \
	    --bam_file sample.bam \
        --reference ~/GRCh37/human_g1k_v37.fa \
        --chromosomes ~/GRCh37/ \
        --samtools_path /opt/tools/samtools/bin/samtools \
        --lumpy_path /home/mcouse/projects/rrg-frid/IMAGINE/SOFTWARE/lumpy-sv/ \
        --lumpy /opt/tools/lumpy/lumpy \
        --CNVnator /opt/tools/CNVnator/src/cnvnator \
        --cnvnator2VCF /opt/tools/CNVnator-0.3.3/cnvnator2VCF.pl \
        --working_dir ~/sample/ \
        --output_prefix sample \
        --read_length 150 \
        --min_mapq 20 \
        --min_support 5 \
        --manta_config /opt/tools/manta/bin/configManta.py \
        --ERDS_path /opt/tools/erds1.1/erds_pipeline.pl \
        --small_variant_vcf  sample.hc.vcf.gz \
        --SV_calling_only

```

### Annotation only

```

python everythingSV.py \
        --bam sample.bam \
        --reference ~/GRCh37/human_g1k_v37.fa \
        --samtools_path /opt/tools/samtools/bin/samtools  \
        --working_dir ~/sample/ \
        --output_prefix sample \
        --SURVIVOR_path /opt/tools/SURVIVOR-src/Debug/SURVIVOR \
        --SURVIVOR_min_support 1 \
        --SURVIVOR_type 1 \
        --SURVIVOR_strand 0 \
        --SURVIVOR_estimate_distance 0 \
        --SURVIVOR_min_SV_small 50 \
        --SURVIVOR_min_SV_large 5000 \
        --SURVIVOR_max_dist_small 200 \
        --SURVIVOR_max_dist_large 1000 \
        --SURVIVOR_additional_vcfs "sample.DELLY.DEL.vcf" \
        "sample.DELLY.DUP.vcf" \
        "sample.DELLY.INV.vcf" \
        --annotation_only \
        --LUMPY_vcf sample/sample.LUMPY.vcf \
         --CNVnator_vcf sample/sample.CNVcall.1000.filter.vcf \
        --Manta_vcf sample/sampleP_Manta/results/variants/diploidSV.vcf \
        --ERDS_vcf sample/sample.ERDS/CP003-P.erds.vcf \
        --table_annovar /opt/tools/annovar/table_annovar.pl \
        --humandb /opt/tools/annovar/humandb \
        --annovar_reference hg19 \
        --deletion_bedfiles "DGV_goldstandard_del.bed, 8, 0.9, TRUE, DGV" \
        "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" \
        --duplication_bedfiles "DGV_goldstandard_dup.bed, 8, 0.9, TRUE, DGV" \
        "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" \
        --inversion_bedfiles "1000G_inv.bed, 8, 0.9, TRUE, 1000G" \
        "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" \
        --insertion_bedfiles "1000G_ins.bed, 8, 0.9, TRUE, 1000G"" \
        "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" \
        "total.combined.domain.named, 4, NA, FALSE, TAD_boundary" \
        --breakend_bedfiles "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" 

```

### SV calling and Annotation

```
python everythingSV.py \
        --bam sample.bam \
        --reference ~/GRCh37/human_g1k_v37.fa \
        --samtools_path /opt/tools/samtools/bin/samtools  \
        --working_dir ~/sample/ \
        --output_prefix sample \
        --lumpy_path /home/mcouse/projects/rrg-frid/IMAGINE/SOFTWARE/lumpy-sv/ \
        --lumpy /opt/tools/lumpy/lumpy \
        --CNVnator /opt/tools/CNVnator/src/cnvnator \
        --cnvnator2VCF /opt/tools/CNVnator-0.3.3/cnvnator2VCF.pl \
        --read_length 150 \
        --min_mapq 20 \
        --min_support 5 \
        --manta_config /opt/tools/manta/bin/configManta.py \
        --ERDS_path /opt/tools/erds1.1/erds_pipeline.pl \
        --small_variant_vcf  sample.hc.vcf.gz \
        --SURVIVOR_path /opt/tools/SURVIVOR-src/Debug/SURVIVOR \
        --SURVIVOR_min_support 1 \
        --SURVIVOR_type 1 \
        --SURVIVOR_strand 0 \
        --SURVIVOR_estimate_distance 0 \
        --SURVIVOR_min_SV_small 50 \
        --SURVIVOR_min_SV_large 5000 \
        --SURVIVOR_max_dist_small 200 \
        --SURVIVOR_max_dist_large 1000 \
        --SURVIVOR_additional_vcfs "sample.DELLY.DEL.vcf" \
        "sample.DELLY.DUP.vcf" \
        "sample.DELLY.INV.vcf" \
        --table_annovar /opt/tools/annovar/table_annovar.pl \
        --humandb /opt/tools/annovar/humandb \
        --annovar_reference hg19 \
        --deletion_bedfiles "DGV_goldstandard_del.bed, 8, 0.9, TRUE, DGV" \
        "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" \
        --duplication_bedfiles "DGV_goldstandard_dup.bed, 8, 0.9, TRUE, DGV" \
        "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" \
        --inversion_bedfiles "1000G_inv.bed, 8, 0.9, TRUE, 1000G" \
        "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" \
        --insertion_bedfiles "1000G_ins.bed, 8, 0.9, TRUE, 1000G"" \
        "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" \
        "total.combined.domain.named, 4, NA, FALSE, TAD_boundary" \
        --breakend_bedfiles "OMIM_genemap2_hg19.bed, 4, NA, FALSE, OMIM" 
```

