# Modified_EcoTyper_Code

### EcoTyper Overview
EcoTyper is a pipeline that allows for the discovery of co-associated cell states based on sample-wide abundance values. For the scRNA-seq pipeline, EcoTyper traditionally performs non-negative matrix factorization (NMF) on a cell type-specific scRNA-seq expression matrix. Using the comphenetic cutoff, EcoTyper identifies an optimal rank k, representing the number of de novo cell states discovered for the given cell type. Similarly, NMF outputs a basis W and coefficient H matrix. W represents the average gene expression per each discovered state and is used for downstream recovery of the k cell states in bulk RNA-seq data. H represents the characteristic abundance of each cell state across all cells of a given cell type; this allows for labeling of the cell state based on the max abundance value for each cell. 

The discovery of de novo cell states potentially inhibits the substitution of labeled cell states. This modification accounts for more granular cell state labeling and recovery. 

### File Setup

The code base is organized based on the specific files used in the modified EcoTyper pipeline; other functional portions of EcoTyper not used specifically in the modified version have been omitted or stored in separate folders. Any file or script used in the modified version contains the suffix, “_Predefined_States_Mode” or “_PDSM.

The scRNA-seq ecotype discovery pipeline has three primary associated files:
  1. EcoTyper_discovery_scRNA_PDSM_Compiled_Steps.R: Calls each step required for ecotype discovery; the subfiles are stored in a separate folder and referenced within the set of compiled steps and detailed below. 
  2. EcoTyper_scRNA_Discovery_PDSM_Config_File.yml: The configuration file takes the file paths of the scRNA-seq expression to be decomposed and its respective metadata file.
  3. EcoTyper_scRNA_Discovery_PDSM_Script.sh: Calls (1) using file path inputs from (2). 
