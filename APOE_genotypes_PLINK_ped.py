#!/bin/env python

# Determine APOE genotypes from PLINK output
    # January 2021
    # Mary B. Makarious and Makayla Portley (LNG/NIA/NINDS/NIH)
    # Script usage:
        # python APOE_genotypes_PLINK_ped.py -i INPUT.ped -o OUTPUT_NAME

# Pre-processing that needs to be done before using this script!
## Pull out the 2 variants of interest
    # plink --bfile YOUR_FILE --snps rs429358,rs7412 --make-bed --out apoe_snps

## Recode into compound genotypes (these will be in the order of how you extracted the SNPs)
    # plink --bfile apoe_snps --recode compound-genotypes --out apoe_snps
    # This makes a .ped and .map file (with no headers!)

    # Format example
    # FID IID PAT MAT SEX PHENO rs429358 rs7412
    # BF-1009 BF-1009 0 0 1 1 CT CC
    # BF-1010 BF-1010 0 0 1 2 TT CC

## APOE Information
# https://www.snpedia.com/index.php/APOE

    # |          APOE GENO         | rs429358 | rs7412 |             COMBINED             |
    # |:--------------------------:|:--------:|:------:|:--------------------------------:|
    # |            e2/e2           |    TT    |   TT   |               TT_TT              |
    # |            e2/e3           |    TT    |   TC   |          TT_TC or TT_CT          |
    # | e2/e4 or e1/e3 (Ambiguous) |    TC    |   TC   | TC_TC or CT_CT or TC_CT or CT_TC |
    # |            e3/e3           |    TT    |   CC   |               TT_CC              |
    # |            e3/e4           |    TC    |   CC   |          TC_CC or CT_CC          |
    # |            e4/e4           |    CC    |   CC   |               CC_CC              |

# Import the necessary packages
import numpy as np
import pandas as pd
import sys
import argparse

# Initialize parser and add arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help="Input file name (with suffix)")
parser.add_argument("--output", "-o", help="Desired output name (without suffix)")
args = parser.parse_args()

# Read in the .ped file and force column names
header_text = ["FID", "IID", "PAT", "MAT", "SEX", "PHENO", "rs429358", "rs7412"]
input_ped_df = pd.read_csv(args.input, sep = " ", header=None, names=header_text)

# Make a combined column, gluing the genotypes from the rs429358 and rs7412 columns
input_ped_df['rs429358_rs7412'] = input_ped_df['rs429358'].astype(str)+'_'+input_ped_df['rs7412']

# Initialize a dictionary with the genotypes to search what genotype the alleles generate
apoe_genotypes_dict = {
    'TT_TT' : 'e2/e2',
    'TT_TC' : 'e2/e3',
    'TT_CT' : 'e2/e3',
    'TC_TC' : 'e2/e4 or e1/e3',
    'CT_CT' : 'e2/e4 or e1/e3',
    'TC_CT' : 'e2/e4 or e1/e3',
    'CT_TC' : 'e2/e4 or e1/e3',
    'TT_CC' : 'e3/e3',
    'TC_CC' : 'e3/e4',
    'CT_CC' : 'e3/e4',
    'CC_CC' : 'e4/e4'
}

# Map the combined column to the dictionary to extract the genotypes
input_ped_df['APOE_GENOTYPE'] = input_ped_df['rs429358_rs7412'].map(apoe_genotypes_dict)

# If any of the combined alleles weren't in the dictionary, the dataframe now has NaN values
# This happens if you have a 0 or missingness somewhere, resulting in an unsure genotype call
# Replace these with something more useful, and state the APOE genotype as "unknown"
input_ped_df.replace(np.nan, 'unknown', regex=True, inplace=True)

# Make a file of just the FID, IID, SEX, PHENO, and APOE genotype
subset_geno_df = input_ped_df.drop(columns=['PAT', 'MAT', 'rs429358', 'rs7412'])

# Print the genotype counts to the user if being used interactively
print("These are the final counts of the number of genotypes")
print(subset_geno_df['APOE_GENOTYPE'].value_counts())

## Export
complete_df_output = args.output + ".APOE_GENOTYPES.csv"
counts_df_output = args.output + ".APOE_SUMMARY.csv"

# Save out the complete dataframe as a .csv
print(f"Your complete genotype file has been saved here: {complete_df_output}")
subset_geno_df.to_csv(complete_df_output, index=False)

# Save out the counts as a .csv
print(f"The summary counts have been saved here: {counts_df_output}")
subset_geno_df['APOE_GENOTYPE'].value_counts().reset_index().to_csv(counts_df_output, index=False, header=['APOE_GENOTYPE', 'COUNT'])

# Done!
print("Thanks!")
