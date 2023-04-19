suppressPackageStartupMessages({
library(config)
library(argparse)
source("pipeline/lib/config.R") 
source("pipeline/lib/misc.R") 
source("pipeline/lib/multithreading.R")
source("pipeline/lib/ecotyper_NMF_Generate_W_Function_Modified.R")
})

parser <- ArgumentParser(add_help = F)

arguments = parser$add_argument_group('Arguments')

arguments$add_argument("-c", "--config", type = "character", metavar="<PATH>", 
    help="Path to the config files [required].")
arguments$add_argument("-h", "--help", action='store_true', help="Print help message.")

args <- parser$parse_args()

if(args$h || is.null(args$config))
{
	parser$print_help()
	quit()
}

config_file = abspath(args$config)

config <- config::get(file = config_file)
# check_discovery_configuration_scRNA(config)

discovery = config$Input$"Discovery dataset name"
discovery_type = config$Input$"Expression type"
scale_column = config$Input$"Annotation file column to scale by"
additional_columns = config$Input$"Annotation file column(s) to plot"

final_output = config$"Output"$"Output folder"

n_threads = config$"Pipeline settings"$"Number of threads"

nmf_restarts = config$"Pipeline settings"$"Number of NMF restarts"
max_clusters = config$"Pipeline settings"$"Maximum number of states per cell type"
cohpenetic_cutoff = config$"Pipeline settings"$"Cophenetic coefficient cutoff"
skip_steps = config$"Pipeline settings"$"Pipeline steps to skip"
p_value_cutoff = config$"Pipeline settings"$"Jaccard matrix p-value cutoff"

# Custom inputs for EcoTyper PDSM.
annotation_file_path = config$"PDSM Settings"$"Annotation File"
expression_matrix_file_path = config$"PDSM Settings"$"Expression Matrix"
input_folder_file_path = config$"PDSM Settings"$"Input Folder"
states_output_folder = config$"PDSM Settings"$"States Output Folder"
ecotypes_output_folder = config$"PDSM Settings"$"Ecotypes Output Folder"
p_val_cutoff = config$"PDSM Settings"$"P-Value Cutoff"
min_states = config$"PDSM Settings"$"Minimum States Per Ecotype"

suppressWarnings({
	final_output = abspath(final_output)	
})

# Starting EcoTyper.
setwd("pipeline")
start = Sys.time()

if(config$"Pipeline settings"$"Filter non cell type specific genes")
{
	fractions = "Cell_type_specific_genes"
}else{
	fractions = "All_genes"
}

# STEP 1: Extract cell type specific genes.

if(!1 %in% skip_steps & config$"Pipeline settings"$"Filter non cell type specific genes")
{
	cat("\nStep 1 (extract cell type specific genes)...\n")

    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
	for(cell_type in cell_types)
	{
		print(cell_type)

        PushToJobQueue(paste("Rscript Pipeline/S1_state_discovery_scRNA_filter_genes_Predefined_States_Mode.R", cell_type, annotation_file_path, expression_matrix_file_path, states_output_folder))	
        
	}
	RunJobQueue()	
	cat("Step 1 (extract cell type specific genes) finished successfully!\n")
	
}else{
	cat("Skipping step 1 (extract cell type specific genes)...\n")
}


# STEP 2: Filtering cell type specific genes. Note, state discovery NMF and the combining restarts steps are skipped.

if(!2 %in% skip_steps)
{
    cat("\nStep 2 (gene filtering on each cell type): filtering cell type specific genes...\n")

    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

	for(cell_type in cell_types)
	{
        PushToJobQueue(paste("Rscript Pipeline/S2_state_discovery_scRNA_distances_Predefined_States_Mode.R", cell_type, TRUE, annotation_file_path, expression_matrix_file_path, input_folder_file_path, states_output_folder))	
        
	}	
	RunJobQueue() 
	
        
}else{
	cat("Skipping step 2 (cell state discovery on correrlation matrices)...\n")
}	

# STEP 5: Create metadata files and coefficient matrix (W) from scRNA-seq data.

# STEP 5 P1: Create metadata files.  

# Define input variables for Python script.
annotationPath = annotation_file_path
sampleColumn = "Sample"
CellTypeColumn = "CellType"
cellStateColumn = "State"
outputPath = states_output_folder

if(!51 %in% skip_steps) # NOTE STEP 5 P1 = 51
{
	# cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")
    cat("\nStep 5 P1 (generating relevant files for ecotype discovery)...\n")

    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

    PushToJobQueue(paste("python Pipeline/S5_P1_EcoTyper_scRNA_Discovery_File_Generation_Predefined_States_Mode.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	

	RunJobQueue() 
	    
}else{
	cat("Skipping step 5 P1 (generating relevant files for ecotype discovery)...\n")
}	


# STEP 5 P2: Create cell state abundances file based on labeled metadata.

if(!52 %in% skip_steps) # NOTE STEP 5 P2 = 52
{
    cat("\nStep 5 P2 (generating state abundances for each cell type)...\n")

      annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))

    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

    PushToJobQueue(paste("python Pipeline/S5_P2_EcoTyper_scRNA_Discovery_Generate_State_Abundances_Predefined_States_Python.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	

	RunJobQueue() 
	    
}else{
	cat("Skipping step 5 P2 (generating state abundances for each cell type)...\n")
}	

# STEP 5 P3: Create binary cell state assignment basis matrix (H) based on labeled metadata for NMF input in step 5 P5. 

if(!53 %in% skip_steps) # NOTE STEP 5 P3 = 53
{
    cat("\nStep 5 P3 (generating binary cell state matrix H)...\n")

    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

    PushToJobQueue(paste("python Pipeline/S5_P3_EcoTyper_scRNA_Discovery_Generate_Binary_H_Predefined_States_Python.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	

	RunJobQueue() 
	    
}else{
	cat("Skipping step 5 P3 (generating binary cell state matrix H)...\n")
}	

# STEP 5 P4: Extract cell state features

if(!54 %in% skip_steps) #NOTE, IT IS STEP 5 P4 = 54
{
	cat("\nStep 5 P4 (extracting cell state features)...\n")
    
    key = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", states_output_folder, "rank_data.txt"))
    
	for(cell_type in key[,1])
	{	
        cat(paste(cell_type, "\n"))
		cat(paste("Extracting marker genes for cell states defined in:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
        
        PushToJobQueue(paste("Rscript Pipeline/S5_P4_state_discovery_extract_features_scRNA_Predefined_States_Mode.R", cell_type, annotation_file_path, expression_matrix_file_path, input_folder_file_path, states_output_folder, n_clusters))
	}	
	RunJobQueue()
	
}else{
	cat("Skipping step 5 P4 (extracting cell state features)...\n")
}	

# STEP 5 P5: Conduct NMF using basis matrix (H) created in STEP 5 P3 to generate coefficient matrix (W). 

if(!55 %in% skip_steps) # NOTE STEP 5 P5 = 55
{
    cat("\nStep 5 P5 (generating average gene expression matrix W)...\n")

    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
    data_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", states_output_folder)
    
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
    
	for(cell_type in cell_types)
	{
        cat(paste0("\nConducting Regular NMF for ", cell_type, "\n"))

        # DOWNSIZING CONTRAINTS FOR LARGE MATRIX - CREATE FUNCTION FOR DOWNSIZING
#         scaledExpMat = read.delim(file.path(data_dir, cell_type, "expression_top_genes_scaled_step_5P1_DS.txt"), row.names = 1)
#         cellTypeBinaryH = read.delim(file.path(data_dir, cell_type, "Binary_H_DS.txt"), row.names = 1)

        scaledExpMat = read.delim(file.path(data_dir, cell_type, "expression_top_genes_scaled_step_5P1.txt"))
        cellTypeBinaryH = read.delim(file.path(data_dir, cell_type, "Binary_H.txt"), row.names = 1)
        avgGeneExpMatW = NMFGenerateW(cellTypeBinaryH, scaledExpMat)

        write.table(avgGeneExpMatW, file.path(data_dir, cell_type, "Average_Cell_State_Gene_Expression_Matrix_W.txt"), sep = "\t")
        
	}	
	RunJobQueue() 
	    
}else{
	cat("Skipping step 5 P5 (generating average gene expression matrix W)...\n")
}	


# STEP 6: Generate gene info files.

if(!6 %in% skip_steps) 
{
    cat("\nStep 6 (generating initial gene info files)...\n")
    
    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
    
	for(cell_type in cell_types)
	{
        
        cat(paste0("\nGenerating initial gene info text files for ", cell_type, "\n"))
                
        PushToJobQueue(paste("Rscript Pipeline/S6_state_discovery_initial_plots_Predefined_States_Mode.R", cell_type, states_output_folder, annotation_file_path))	
        
	}	
	RunJobQueue() 
	    
}else{
	cat("Skipping Step 6 (generating initial gene info files)...\n")
}	

# STEP 8: Ecotype discovery.

# STEP 8 P1: Identify ecotypes.

if(!81 %in% skip_steps)
{
    cat("\nStep 8 P1 (ecotype discovery)...\n")
    
    PushToJobQueue(paste("Rscript Pipeline/S8_P1_ecotypes_scRNA_PDSM_Updated.R", states_output_folder, annotation_file_path, p_val_cutoff, min_states)) 
    
	RunJobQueue()
    
    cat("Step 8 P1 (ecotype discovery) finished successfully!\n")

}else{
	cat("\nSkipping step 8 P1 (ecotype discovery)...\n")
}

# STEP 8 P2: Assign ecotypes.

if(!82 %in% skip_steps)
{
    cat("Step 8 P2 (ecotype assignment) assigning ecotypes...!\n")
    
    PushToJobQueue(paste("Rscript Pipeline/S8_P2_ecotypes_assign_samples_scRNA_PDSM_Updated.R", states_output_folder, "State", annotation_file_path, paste(additional_columns, collapse = " "))) 
	
	RunJobQueue()
    cat("Step 8 P2 (ecotype assignment) ecotype assignment finished successfully!\n")

}else{
	cat("\nSkipping step 8 P2 (ecotype assignment)...\n")
}


# cat("\nCopying EcoTyper results to the output folder!\n")

# if(file.exists(final_output) && length(list.files(final_output)) > 0)
# {
# 	old_results_folder = paste0(final_output, format(Sys.time(), " %a %b %d %X %Y"))
# 	dir.create(old_results_folder, recursive = T, showWarnings = F)
# 	warning(paste0("The output folder contains files from a previous run. Moving those files to: '", old_results_folder, "'"))	
# 	system(paste0("mv -f ", final_output, "/* '", old_results_folder, "'"))
# }

# dir.create(final_output, recursive = T, showWarnings = F)

# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"), final_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_plot.pdf"), final_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_plot.png"), final_output))

# key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# for(cell_type in key[,1])
# {		
# 	n_clusters = key[key[,1] == cell_type, 2]	
# 	ct_output = file.path(final_output, cell_type)
# 	dir.create(ct_output, recursive = T, showWarnings = F)
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "gene_info.txt"), ct_output))
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_abundances.txt"), ct_output))
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment.txt"), ct_output))
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment_heatmap.pdf"), ct_output))
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "state_assignment_heatmap.png"), ct_output))
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "heatmap_data.txt"), ct_output))
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "heatmap_top_ann.txt"), ct_output))	
# }	

# ct_output = file.path(final_output, "Ecotypes")
# dir.create(ct_output, recursive = T, showWarnings = F)
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotypes.txt"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotype_assignment.txt"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "ecotype_abundance.txt"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "heatmap_assigned_samples_viridis.pdf"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "heatmap_assigned_samples_viridis.png"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "jaccard_matrix.pdf"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "jaccard_matrix.png"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "nclusters_jaccard.png"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "discovery", "nclusters_jaccard.pdf"), ct_output))

# Ending EcoTyper
end = Sys.time()
# cat(paste0("\nEcoTyper finished succesfully! Please find the results in: '", final_output, "'.\nRun time: ", format(end - start, digits = 0), "\n"))

