# Input data:

All the necessary files are in the folder Input_data/
1. The expression data in different files (one file for each condition) in the folder 4565.zip
2. The homeologs list homoeolog_CSv2_ABD.list
3. The files CSv2.1_HC_vs_NCBI_besthit_pair.txt and CSv2.1_LC_vs_NCBI_besthit_pair.txt to go from TraesID to RefSeqID
4. The file GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_cds_from_genomic.fna to go from EGID to RefSeqID and get genes length and isoform number
5. The file iwgsc_refseqv2.1_functional_annotation.csv for genes functional annotation
6. The file run-study-attrib-biosample.txt to have information about the studies the expression data come from. 

# Scripts order:

## Expression Data Preprocessing:
1. Preprocessing/Expression_Data/Merge_expression_data.R
2. Preprocessing/Expression_Data/TPM_UQ_Normalization.R
3. Preprocessing/Expression_Data/Batch_Normalization_ComBat.R
4. Preprocessing/Expression_Data/Filter_genes_zero_var.R
5. Preprocessing/Expression_Data/Visualization_conditions.R

## Homeologs List Preprocessing:
1. To create RefseqID_2_EGI:
perl -lane 'print "$1\t$2" if /^>/ && /(NC_\S+)_.*GeneID:(\d+)/' GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_cds_from_genomic.fna 
2. Preprocessing/Homeologs_list/Create_Traes_2_EGI.py
3. Preprocessing/Homeologs_list/Homeologs_traes_2_EGI.R

## Get genes length and isoform number:
1. Preprocessing/Genes_length_isoform_number.py

## Comparisons and categorization of triads on different properties:
With xx the property: Expression_levels, Expression_profiles, GO_annotations, Length, Isoform_number, Chromosomes
1. xx_comparisons/xx_comparisons.R
2. xx_comparisons/xx_category.py

## Visualization:
1. Visualization/Pie_Plot_script.R

## Crossing the results:
1. Crossing_results/Cross_strategy.R
