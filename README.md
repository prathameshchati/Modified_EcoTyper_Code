# Modified_EcoTyper_Code

## EcoTyper Overview

[EcoTyper](https://github.com/digitalcytometry/ecotyper) is a pipeline that allows for the discovery of co-associated cell states based on sample-wide abundance values. For the scRNA-seq discovery pipeline, EcoTyper traditionally performs non-negative matrix factorization (NMF) on a cell type-specific scRNA-seq expression matrix. Using the comphenetic cutoff, EcoTyper identifies an optimal rank __k__, representing the number of _de novo_ cell states discovered for the given cell type. Similarly, NMF outputs a basis __W__ and coefficient __H__ matrix. __W__ represents the average gene expression per each discovered state and is used for downstream recovery of the __k__ cell states in bulk RNA-seq data. __H__ represents the characteristic abundance of each cell state across all cells of a given cell type; this allows for labeling of the cell state based on the max abundance value for each cell. 

The discovery of de novo cell states potentially inhibits the substitution of labeled cell states. This modification accounts for more granular cell state labeling and recovery. 

## Directory Structure

The code base is organized based on the specific files used in the modified EcoTyper pipeline; other functional portions of EcoTyper not used specifically in the modified version have been omitted or stored in separate folders. Any file or script used in the modified version contains the suffix, "_Predefined_States_Mode" or "_PDSM".
  
### scRNA Discovery

The scRNA discovery code can be found in __scRNA_Discovery_PDSM__. The scRNA-seq ecotype discovery pipeline has three primary associated files:

  1. __EcoTyper_discovery_scRNA_PDSM_Compiled_Steps.R__: Calls each step required for ecotype discovery; the subfiles are stored in a separate folder and referenced within the set of compiled steps and detailed below. 
  2. __EcoTyper_scRNA_Discovery_PDSM_Config_File.yml__: The configuration file takes the file paths of the scRNA-seq expression to be decomposed and its respective metadata file.
  3. __EcoTyper_scRNA_Discovery_PDSM_Script.sh__: Calls __(1)__ using file path inputs from __(2)__. 
  
The subfolder __Pipeline__ contains labeled steps that are called within __(1)__ and beging with __S__ followed by its respective order. The __P__ indicates subparts of a certain step. The steps are not numbered in equal steps since each step corresponds to its true EcoTyper step. Certain steps were excluded based on their necessity within and compatibility with the modified pipeline. Below, each step is described, including its function and respective input and output files.

#### scRNA Discovery Steps

##### S1

