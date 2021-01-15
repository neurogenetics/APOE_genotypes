# !/bin/env bash

# Determine APOE genotypes from PLINK output
    # January 2021
    # Mary B. Makarious, Makayla Portley, and Cornelis Blauwendraat (LNG/NIA/NINDS/NIH)

## On Biowulf
# Load an interactive node

# Load the necessary modules
module load plink #v1.9.0
module load python #3.7

# Initialize workspace variables
WORK_DIR="/path/to/data"

## APOE Information
    # |          APOE GENO         	| rs429358 	| rs7412 	|             COMBINED             	|
    # |:--------------------------:	|:--------:	|:------:	|:--------------------------------:	|
    # |            e1/e1           	|    CC    	|   TT   	|               CC_TT              	|
    # |            e1/e2           	|    CT    	|   TT   	|          CT_TT or TC_TT          	|
    # |            e1/e4           	|    CC    	|   CT   	|          CC_CT or CC_TC          	|
    # |            e2/e2           	|    TT    	|   TT   	|               TT_TT              	|
    # |            e2/e3           	|    TT    	|   TC   	|          TT_TC or TT_CT          	|
    # | e2/e4 or e1/e3 (Ambiguous) 	|    TC    	|   TC   	| TC_TC or CT_CT or TC_CT or CT_TC 	|
    # |            e3/e3           	|    TT    	|   CC   	|               TT_CC              	|
    # |            e3/e4           	|    TC    	|   CC   	|          TC_CC or CT_CC          	|
    # |            e4/e4           	|    CC    	|   CC   	|               CC_CC              	|

cd $WORK_DIR

## Pull out the 2 variants of interest
plink --bfile ./../YOUR_PLINK_FILE --snps rs429358,rs7412 --make-bed --out apoe_snps

## Recode into compound genotypes (these will be in the order of how you extracted the SNPs)
plink --bfile apoe_snps --recode compound-genotypes --out apoe_snps
    # This makes a .ped and .map file (with no headers!)
    # Format example
    # FID IID PAT MAT SEX PHENO rs429358 rs7412
    # sample1 sample1 0 0 1 1 CT CC
    # sample2 sample2 0 0 1 2 TT CC

## Run the following Python script
python APOE_genotypes_PLINK_ped.py -i apoe_snps.ped -o apoe_snps_test
