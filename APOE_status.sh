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
WORK_DIR="/path/to/data"

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
