# Modified_EcoTyper_Code

## EcoTyper Overview

[EcoTyper](https://github.com/digitalcytometry/ecotyper) is a pipeline that allows for the discovery of co-associated cell states based on sample-wide abundance values. For the scRNA-seq discovery pipeline, EcoTyper traditionally performs non-negative matrix factorization (NMF) on a cell type-specific scRNA-seq expression matrix. Using the comphenetic cutoff, EcoTyper identifies an optimal rank __k__, representing the number of _de novo_ cell states discovered for the given cell type. Similarly, NMF outputs a basis __W__ and coefficient __H__ matrix. __W__ represents the average gene expression per each discovered state and is used for downstream recovery of the __k__ cell states in bulk RNA-seq data. __H__ represents the characteristic abundance of each cell state across all cells of a given cell type; this allows for labeling of the cell state based on the max abundance value for each cell. 

The discovery of de novo cell states potentially inhibits the substitution of labeled cell states. This modification accounts for more granular cell state labeling and recovery. 

## Directory Structure

The code base is organized based on the specific files used in the modified EcoTyper pipeline; other functional portions of EcoTyper not used specifically in the modified version have been omitted or stored in separate folders. Any file or script used in the modified version contains the suffix, "_Predefined_States_Mode" or "_PDSM".

The pipeline is split into two parts which are described briefly here and more in detail below:
  1. scRNA Discovery: Identifies ecotypes from labeled scRNA-seq data and generates custom cell state signature matrix for bulk RNA-seq recovery.
  2. Bulk Recovery: Recovers identified cell states and ecotypes in bulk RNA-seq data. 
  
### scRNA Discovery

The scRNA discovery code can be found in __scRNA_Discovery_PDSM__. The scRNA-seq ecotype discovery pipeline has three primary associated files:

  1. __EcoTyper_discovery_scRNA_PDSM_Compiled_Steps.R__: Calls each step required for ecotype discovery; the subfiles are stored in a separate folder and referenced within the set of compiled steps and detailed below. 
  2. __EcoTyper_scRNA_Discovery_PDSM_Config_File.yml__: The configuration file takes the file paths of the scRNA-seq expression to be decomposed and its respective metadata file.
  3. __EcoTyper_scRNA_Discovery_PDSM_Script.sh__: Calls __(1)__ using file path inputs from __(2)__. 
  
The subfolder __Pipeline__ contains labeled steps that are called within __(1)__ and beging with __S__ followed by its respective order. The __P__ indicates subparts of a certain step. The steps are not numbered in equal steps since each step corresponds to its true EcoTyper step. Certain steps were excluded based on their necessity within and compatibility with the modified pipeline. Below, each step is described, including its function and respective input and output files. Note, specific file outputs are specified in the [EcoTyper](https://github.com/digitalcytometry/ecotyper) documentation while file specifically modified here are detailed below. 

#### Steps

__Note:__ 

##### S1: Isolate cell type-specific genes

- File Path: Pipeline/S1_state_discovery_scRNA_filter_genes_Predefined_States_Mode.R
- Function: Identifies top differentially expressed genes per cell type. 
- Inputs: Annotations, expression atrix
- Outputs: Cell type-specific genes (_cell_type_specific_genes_raw.txt, _cell_type_specific_genes.txt)

##### S2: Scales and filters expression matrix

- File Path: Pipeline/S2_state_discovery_scRNA_distances_Predefined_States_Mode.R
- Function: Scales, normalizes, and filters gene expression matrix per cell type; generates cell type-specific gene expression matrices. 
- Inputs: S1 outputs, annotations, expression, matrix
- Outputs: Cell type-specific expression matrices (expression_full_matrix_log2.txt, expression_full_matrix_scaled.txt, expression_top_genes_scaled_filt.txt)

##### S5: File generation and NMF

###### P1: Generates metadata for NMF

- File Path: Pipeline/S5_P1_EcoTyper_scRNA_Discovery_File_Generation_Predefined_States_Mode.py
- Function: First generates mappings between initial and final state labels (the difference in the modified pipeline is arbitrary and can also be the same labels). It next generates the initial and state assignment files for ecotype discovery. Finally, the rank data is generated, which indicates the number of cell states per cell type. The rank should be fixed given each cell type is manually subclustered to a finite number of cell states. 
- Inputs: Annotations (with appropriate cell type and state columns containing labels)
- Outputs: mapping_to_initial_states.txt, initial_state_assignment.txt, state_assignment.txt, rank_data.txt

###### P2: Generates state abundances

- File Path: Pipeline/S5_P2_EcoTyper_scRNA_Discovery_Generate_State_Abundances_Predefined_States_Python.py
- Function: Uses cell state labels and generates fractional abundance of each cell state per cell type across all samples. Note, cell state fractions are scaled to 1 within their respective cell types. 
- Inputs: Annotations (with appropriate cell type and state columns containing labels)
- Outputs: state_abundances.txt

###### P3: Generate binary H matrix

- File Path: Pipeline/S5_P3_EcoTyper_scRNA_Discovery_Generate_Binary_H_Predefined_States_Python.py
- Function: The binary H matrix is fitted to the cell type-specific NMF model to recover W, which gives the average gene expression per labeled cell state. Within each cell type, the cell state with the max abundance is assigned a 1 while the remaining cell states are assigned a 0; this is repeated across all samples. The bimary H matrix is a substitution for the auto-generated H matrix in the traditional EcoTyper pipeline. 
- Inputs: S5 P2 outputs, annotations
- Outputs: Binary_H.txt

###### P4: Retrieves features associated with cell states

- File Path: Pipeline/S5_P4_state_discovery_extract_features_scRNA_Predefined_States_Mode.r
- Function: Computes cell state-specific differentially expressed genes and isolates the expression matrix for each cell type based on the top signature genes for each cell state. This supports the NMF step by limiting the space of possible signature genes to converge on. 
- Inputs: S2 outputs, S5 P1 outputs
- Outputs: full_FCs.txt, top_FCs.txt, expression_top_genes_log2_step_5P1.txt, expression_top_genes_scaled_step_5P1.txt

###### P5: Conduct modified NMF on cell type-specific scRNA-seq matrix

- File Path: Directly within scRNA_Discovery_PDSM/EcoTyper_discovery_scRNA_PDSM_Compiled_Steps.R
- Function: Uses the binary H matrix and fits an NMF model that enables the recovery of specific W which represents the average gene expression per each labeled cell state for a given cell type. 
- Inputs: S5 P3 outputs (Binary_H.txt), S5 P4 outputs (expression_top_genes_scaled_step_5P1.txt)
- Outputs: Average_Cell_State_Gene_Expression_Matrix_W.txt

##### S6: Create gene information file

- File Path: Pipeline/S6_state_discovery_initial_plots_Predefined_States_Mode.R 
- Function: Generates gene information file which contains cell state-specific genes filtered for specific signature genes. 
- Inputs: Annotations, S5 P1 outputs (initial_state_assignment.txt), S5 P4 outputs (expression_top_genes_scaled_step_5P1.txt)
- Outputs: initial_gene_info.txt

##### S8: Ecotype discvoery

###### P1: Identifies ecotypes

- File Path: Pipeline/S8_P1_ecotypes_scRNA_PDSM_Updated.R
- Function: Cell state abundances per cell type are characterized and binarized across all samples. Cell states with similar abundance profiles are clustered to generate communities of co-occuring (co-associated) cell states, termed ecotypes. 
- Inputs: rank_data.txt, mapping_to_initial_states.txt, initial_state_assignment.txt, state_assignment.txt
- Outputs: binary_classification_all_states.txt, jaccard_matrix.txt, silhouette_initial.txt, initial_jaccard_matrix.pdf, ecotypes.txt, silhouette.txt, jaccard_matrix.png

###### P2: Assigns ecotypes

- File Path: Pipeline/S8_P2_ecotypes_assign_samples_scRNA_PDSM_Updated.R
- Function: After ecotypes are discovered, each sample is assigned an ecotype based on the ecotype abundance data. Ecotype abundance heatmap and clustermaps are also generated.
- Inputs: rank_data.txt, ecotypes.txt, mapping_all_states.txt, mapping_to_initial_states.txt, initial_state_assignment.txt
- Outputs:combined_state_abundances.txt, ecotype_abundance.txt, assignment_p_vals.txt, heatmap_assigned_samples_viridis.pdf, heatmap_assigned_samples_viridis.png, heatmap_assigned_samples_YlGnBu.pdf, heatmap_assigned_samples_YlGnBu.png


