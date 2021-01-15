# !/bin/env bash

# Determine APOE genotypes from PLINK output
    # January 2021
    # Mary B. Makarious and Makayla Portley (LNG/NIA/NINDS/NIH)

## On Biowulf
# Load an interactive node
sinteractive --cpus-per-task=2 --mem=200g --constraint=ibfdr --gres=lscratch:200 --time=24:00:00

# Load the necessary modules
module load plink #v1.9.0
module load python

# Initialize workspace variables
WORK_DIR="/data/LNG/makariousmb/Projects/GenoMLxAMP/clean_wgs/exp_apoe"

## APOE Information
    # |  APOE GENO  | rs429358 | rs7412 |             Combined             |
    # |:-----------:|:--------:|:------:|:--------------------------------:|
    # |    e2/e2    |    TT    |   TT   |               TT_TT              |
    # |    e2/e3    |    TT    |   TC   |          TT_TC or TT_CT          |
    # |    e2/e4    |    TC    |   TC   | TC_TC or CT_CT or TC_CT or CT_TC |
    # |    e3/e3    |    TT    |   CC   |               TT_CC              |
    # |    e3/e4    |    TC    |   CC   |          TC_CC or CT_CC          |
    # |    e4/e4    |    CC    |   CC   |               CC_CC              |

cd $WORK_DIR

## Pull out the 2 variants of interest
plink --bfile ./../AMP_Euro_sampleQC_variantQC_FINAL --snps rs429358,rs7412 --make-bed --out apoe_snps

## Recode into compound genotypes (these will be in the order of how you extracted the SNPs)
plink --bfile apoe_snps --recode compound-genotypes --out apoe_snps
    # This makes a .ped and .map file (with no headers!)
    # Format example
    # FID IID PAT MAT SEX PHENO rs429358 rs7412
    # BF-1009 BF-1009 0 0 1 1 CT CC
    # BF-1010 BF-1010 0 0 1 2 TT CC
    # BF-1011 BF-1011 0 0 1 2 CT TC
    # BF-1012 BF-1012 0 0 2 1 TT TC
    # BF-1013 BF-1013 0 0 1 1 TT CC
    # BF-1014 BF-1014 0 0 2 1 CT CC
    # BF-1015 BF-1015 0 0 1 2 TT TC

## Run the following Python script
python APOE_genotypes_PLINK_ped.py -i apoe_snps.ped -o apoe_snps_test
