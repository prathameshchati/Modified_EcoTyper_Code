suppressPackageStartupMessages({
library(config)
library(argparse)
source("pipeline/lib/config.R") 
source("pipeline/lib/misc.R") 
source("pipeline/lib/multithreading.R")
# source("pipeline/lib/ecotyper_NMF_Generate_W_Function.R")
# source("pipeline/lib/ecotyper_NMF_Generate_W_Function_One_State.R")
source("pipeline/lib/ecotyper_NMF_Generate_W_Function_Modified.R")
})

parser <- ArgumentParser(add_help = F)

arguments = parser$add_argument_group('Arguments')

arguments$add_argument("-c", "--config", type = "character", metavar="<PATH>", 
    help="Path to the config files [required].")
arguments$add_argument("-h", "--help", action='store_true', help="Print help message.")

args <- parser$parse_args()
#print(args)
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

# CUSTOM INPUTS
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

#Starting EcoTyper
setwd("pipeline")
start = Sys.time()

if(config$"Pipeline settings"$"Filter non cell type specific genes")
{
	fractions = "Cell_type_specific_genes"
}else{
	fractions = "All_genes"
}

# STEP 1

if(!1 %in% skip_steps & config$"Pipeline settings"$"Filter non cell type specific genes")
{
	cat("\nStep 1 (extract cell type specific genes)...\n")

	# annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
#     annotation = read.delim(file.path(file.path("../example_data", "scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))	
    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
	for(cell_type in cell_types)
	{
		print(cell_type)
		# PushToJobQueue(paste("Rscript state_discovery_scRNA_filter_genes.R", discovery, fractions, cell_type, scale_column))	
#         PushToJobQueue(paste("Rscript state_discovery_scRNA_filter_genes_Predefined_States_Mode.R", cell_type))	
        PushToJobQueue(paste("Rscript state_discovery_scRNA_filter_genes_Predefined_States_Mode.R", cell_type, annotation_file_path, expression_matrix_file_path, states_output_folder))	
        
	}
	RunJobQueue()	
	cat("Step 1 (extract cell type specific genes) finished successfully!\n")
	
}else{
	cat("Skipping step 1 (extract cell type specific genes)...\n")
}


# if(!1 %in% skip_steps & config$"Pipeline settings"$"Filter non cell type specific genes")
# {
# 	cat("\nStep 1 (extract cell type specific genes)...\n")

# 	annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
# 	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
# 	for(cell_type in cell_types)
# 	{
# 		print(cell_type)
# 		PushToJobQueue(paste("Rscript state_discovery_scRNA_filter_genes.R", discovery, fractions, cell_type, scale_column))	
# 	}
# 	RunJobQueue()	
# 	cat("Step 1 (extract cell type specific genes) finished successfully!\n")
	
# }else{
# 	cat("Skipping step 1 (extract cell type specific genes)...\n")
# }

# STEP 2

if(!2 %in% skip_steps)
{
	# cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")
    cat("\nStep 2 (gene filtering on each cell type): filtering cell type specific genes...\n")

	# annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
#     annotation = read.delim(file.path(file.path("../example_data", "scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))	
    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

	for(cell_type in cell_types)
	{
		#filter_genes = (fractions == "Cell_type_specific_genes")
		# PushToJobQueue(paste("Rscript state_discovery_scRNA_distances.R", discovery, fractions, cell_type, filter_genes, scale_column))	
        PushToJobQueue(paste("Rscript state_discovery_scRNA_distances_Predefined_States_Mode.R", cell_type, TRUE, annotation_file_path, expression_matrix_file_path, input_folder_file_path, states_output_folder))	
        
	}	
	RunJobQueue() 
	
    
    if(!25 %in% skip_steps)
    {
        cat("Step 2 (cell state discovery on correrlation matrices): Running NMF (Warning: This step might take a long time!)...\n")

        for(cell_type in cell_types)
        {		
            if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, "expression_top_genes_scaled.txt")))
            {			
                next
            }
            for(n_clusters in 2:max_clusters)
            {
                for(restart in 1:nmf_restarts)
                {
                    if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, n_clusters, "restarts", restart, "estim.RData")))
                    {
                        PushToJobQueue(paste("Rscript state_discovery_NMF.R", "discovery_cross_cor", discovery, fractions, cell_type, n_clusters, restart))
                    }else{					
                        cat(paste0("Warning: Skipping NMF on '", cell_type, "' (number of states = ", n_clusters, ", restart ", restart, "), as the output file '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "restarts", restart, "estim.RData"), "' already exists!\n"))
                    }
                } 
            }			
        } 
        RunJobQueue()

        cat("Step 2 (cell state discovery on correrlation matrices): Aggregating NMF results...\n")
        for(cell_type in cell_types)
        {	
            if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, "expression_top_genes_scaled.txt")))
            {			
                next
            }				
            PushToJobQueue(paste("Rscript state_discovery_combine_NMF_restarts.R", "discovery_cross_cor", discovery, fractions, cell_type, max_clusters, nmf_restarts))
        } 
        RunJobQueue()
        cat("Step 2 (cell state discovery on correrlation matrices) finished successfully!\n")
    }
    
}else{
	cat("Skipping step 2 (cell state discovery on correrlation matrices)...\n")
}	


# if(!2 %in% skip_steps)
# {
# 	cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")

# 	annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
# 	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

# 	for(cell_type in cell_types)
# 	{
# 		filter_genes = (fractions == "Cell_type_specific_genes")
# 		PushToJobQueue(paste("Rscript state_discovery_scRNA_distances.R", discovery, fractions, cell_type, filter_genes, scale_column))	
# 	}	
# 	RunJobQueue() 
	
# 	cat("Step 2 (cell state discovery on correrlation matrices): Running NMF (Warning: This step might take a long time!)...\n")

# 	for(cell_type in cell_types)
# 	{		
# 		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, "expression_top_genes_scaled.txt")))
# 		{			
# 			next
# 		}
# 		for(n_clusters in 2:max_clusters)
# 		{
# 			for(restart in 1:nmf_restarts)
# 			{
# 				if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, n_clusters, "restarts", restart, "estim.RData")))
# 				{
# 					PushToJobQueue(paste("Rscript state_discovery_NMF.R", "discovery_cross_cor", discovery, fractions, cell_type, n_clusters, restart))
# 				}else{					
# 					cat(paste0("Warning: Skipping NMF on '", cell_type, "' (number of states = ", n_clusters, ", restart ", restart, "), as the output file '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "restarts", restart, "estim.RData"), "' already exists!\n"))
# 				}
# 			} 
# 		}			
# 	} 
# 	RunJobQueue()
		
# 	cat("Step 2 (cell state discovery on correrlation matrices): Aggregating NMF results...\n")
# 	for(cell_type in cell_types)
# 	{	
# 		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery_cross_cor", cell_type, "expression_top_genes_scaled.txt")))
# 		{			
# 			next
# 		}				
# 		PushToJobQueue(paste("Rscript state_discovery_combine_NMF_restarts.R", "discovery_cross_cor", discovery, fractions, cell_type, max_clusters, nmf_restarts))
# 	} 
# 	RunJobQueue()
# 	cat("Step 2 (cell state discovery on correrlation matrices) finished successfully!\n")
# }else{
# 	cat("Skipping step 2 (cell state discovery on correrlation matrices)...\n")
# }	

# STEP 3

# if(!3 %in% skip_steps)
# {
# 	cat("\nStep 3 (choosing the number of cell states)...\n")
# 	PushToJobQueue(paste("Rscript state_discovery_rank_selection.R", "discovery_cross_cor", discovery, fractions, max_clusters, cohpenetic_cutoff))
# 	RunJobQueue()
# 	cat("Step 3 (choosing the number of cell states) finished successfully!\n")
# }else{
# 	cat("Skipping step 3 (choosing the number of cell states)...\n")
# }

# if(!3 %in% skip_steps)
# {
# 	cat("\nStep 3 (choosing the number of cell states)...\n")
# 	PushToJobQueue(paste("Rscript state_discovery_rank_selection.R", "discovery_cross_cor", discovery, fractions, max_clusters, cohpenetic_cutoff))
# 	RunJobQueue()
# 	cat("Step 3 (choosing the number of cell states) finished successfully!\n")
# }else{
# 	cat("Skipping step 3 (choosing the number of cell states)...\n")
# }

# STEP 4

# if(!4 %in% skip_steps)
# {
# 	cat("\nStep 4 (extracting cell state information)...\n")

# 	system(paste("cp -f ", config_file, file.path("../EcoTyper", discovery, "config_used.yml")))

# 	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# 	for(cell_type in key[,1])
# 	{	
# 		cat(paste("Extracting cell states information for:", cell_type, "\n"))
# 		n_clusters = key[key[,1] == cell_type, 2]
# 		PushToJobQueue(paste("Rscript state_discovery_initial_plots_scRNA.R", "discovery_cross_cor", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
# 	}	
# 	RunJobQueue()
# 	cat("Step 4 (extracting cell state information) finished successfully!\n")
# }else{
# 	cat("\nSkipping step 4 (extracting cell state information)...\n")
# }

# if(!4 %in% skip_steps)
# {
# 	cat("\nStep 4 (extracting cell state information)...\n")

# 	system(paste("cp -f ", config_file, file.path("../EcoTyper", discovery, "config_used.yml")))

# 	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# 	for(cell_type in key[,1])
# 	{	
# 		cat(paste("Extracting cell states information for:", cell_type, "\n"))
# 		n_clusters = key[key[,1] == cell_type, 2]
# 		PushToJobQueue(paste("Rscript state_discovery_initial_plots_scRNA.R", "discovery_cross_cor", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
# 	}	
# 	RunJobQueue()
# 	cat("Step 4 (extracting cell state information) finished successfully!\n")
# }else{
# 	cat("\nSkipping step 4 (extracting cell state information)...\n")
# }

# STEP 5

# P1

# DEFINE INPUT VARIABLES FOR PYTHON FUNCTION

# anotationPath = "../example_data/scRNA_CRC_Annotation_Ecotype_States_Filtered.txt"
annotationPath = annotation_file_path

# dataPath = "../example_data/scRNA_CRC_data_Filtered.txt"
sampleColumn = "Sample"
CellTypeColumn = "CellType"
cellStateColumn = "State"

# outputPath = "/duo4/users/pchati/ecotyper/Generate_Cell_State_Abundances_Predefined_States_Test_Outputs/Test_5_Example_scRNA_seq_Data"
outputPath = states_output_folder

if(!51 %in% skip_steps) # NOTE STEP 5 P1 = 51
{
	# cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")
    cat("\nStep 5 P1 (generating relevant files for ecotype discovery)...\n")

	# annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
#     annotation = read.delim(file.path(file.path("../example_data", "scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))	
    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

# 	for(cell_type in cell_types)
# 	{
#         PushToJobQueue(paste("python lib/Python_Scripts/EcoTyper_scRNA_Discovery_File_Generation_Predefined_States_Mode.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	
#         # PushToJobQueue(paste("python lib/Python_Scripts/Test_Script.py", "arg1", "arg2", "arg3"))	
        
# 	}	
    
    #### REMOVED UNNECESSARY LOOP
    PushToJobQueue(paste("python lib/Python_Scripts/EcoTyper_scRNA_Discovery_File_Generation_Predefined_States_Mode.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	

	RunJobQueue() 
	    
}else{
	cat("Skipping step 5 P1 (generating relevant files for ecotype discovery)...\n")
}	


# P2

# DEFINE INPUT VARIABLES FOR PYTHON FUNCTION

# anotationPath = "../example_data/scRNA_CRC_Annotation_Ecotype_States_Filtered.txt"
# annotationPath = annotation_file_path

# dataPath = "../example_data/scRNA_CRC_data_Filtered.txt"
# sampleColumn = "Sample"
# CellTypeColumn = "CellType"
# cellStateColumn = "State"

# outputPath = "/duo4/users/pchati/ecotyper/Generate_Cell_State_Abundances_Predefined_States_Test_Outputs/Test_5_Example_scRNA_seq_Data"
# outputPath = states_output_folder

if(!52 %in% skip_steps) # NOTE STEP 5 P2 = 52
{
	# cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")
    cat("\nStep 5 P2 (generating state abundances for each cell type)...\n")

	# annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
#     annotation = read.delim(file.path(file.path("../example_data", "scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))	
      annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))

    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

# 	for(cell_type in cell_types)
# 	{
#         PushToJobQueue(paste("python lib/Python_Scripts/EcoTyper_scRNA_Discovery_Generate_State_Abundances_Predefined_States_Python.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	
#         # PushToJobQueue(paste("python lib/Python_Scripts/Test_Script.py", "arg1", "arg2", "arg3"))	
        
# 	}	
    
    #### REMOVED UNNECESSARY LOOP
    PushToJobQueue(paste("python lib/Python_Scripts/EcoTyper_scRNA_Discovery_Generate_State_Abundances_Predefined_States_Python.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	

	RunJobQueue() 
	    
}else{
	cat("Skipping step 5 P2 (generating state abundances for each cell type)...\n")
}	

# P3

# DEFINE INPUT VARIABLES FOR PYTHON FUNCTION
# anotationPath = "../example_data/scRNA_CRC_Annotation_Ecotype_States_Filtered.txt"
# dataPath = "../example_data/scRNA_CRC_data_Filtered.txt"
# sampleColumn = "Sample"
# CellTypeColumn = "CellType"
# cellStateColumn = "State"
# outputPath = "/duo4/users/pchati/ecotyper/Generate_Cell_State_Abundances_Predefined_States_Test_Outputs/Test_5_Example_scRNA_seq_Data"

if(!53 %in% skip_steps) # NOTE STEP 5 P3 = 53
{
	# cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")
    cat("\nStep 5 P3 (generating binary cell state matrix H)...\n")

	# annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
#     annotation = read.delim(file.path(file.path("../example_data", "scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))	
    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))

    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	

# 	for(cell_type in cell_types)
# 	{
#         PushToJobQueue(paste("python lib/Python_Scripts/EcoTyper_scRNA_Discovery_Generate_Binary_H_Predefined_States_Python.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	
#         # PushToJobQueue(paste("python lib/Python_Scripts/Test_Script.py", "arg1", "arg2", "arg3"))	
        
# 	}	
    
    
    #### REMOVED UNNECESSARY LOOP
    PushToJobQueue(paste("python lib/Python_Scripts/EcoTyper_scRNA_Discovery_Generate_Binary_H_Predefined_States_Python.py", annotationPath, sampleColumn, CellTypeColumn, cellStateColumn, outputPath))	

	RunJobQueue() 
	    
}else{
	cat("Skipping step 5 P3 (generating binary cell state matrix H)...\n")
}	

# P4

if(!54 %in% skip_steps) #NOTE, IT IS STEP 5 P4 = 54
{
	cat("\nStep 5 P4 (extracting cell state features)...\n")

	# key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
    
#     key = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", "Test_5_Example_scRNA_seq_Data", "rank_data.txt"))
    
    key = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", states_output_folder, "rank_data.txt"))
    
#     cell_types = c('Malignant')
#     cell_types = c('Bcell', 'Malignant')

	for(cell_type in key[,1])
#     for(cell_type in cell_types)
	{	
        cat(paste(cell_type, "\n"))
		cat(paste("Extracting marker genes for cell states defined in:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		# PushToJobQueue(paste("Rscript state_discovery_extract_features_scRNA.R", discovery, fractions, cell_type, n_clusters)) 		 
        PushToJobQueue(paste("Rscript state_discovery_extract_features_scRNA_Predefined_States_Mode.R", cell_type, annotation_file_path, expression_matrix_file_path, input_folder_file_path, states_output_folder, n_clusters))
	}	
	RunJobQueue()
	
}else{
	cat("Skipping step 5 P4 (extracting cell state features)...\n")
}	

# P5

# DEFINE INPUT VARIABLES FOR GENERATE W FUNCTION

if(!55 %in% skip_steps) # NOTE STEP 5 P5 = 55
{
	# cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")
    cat("\nStep 5 P5 (generating average gene expression matrix W)...\n")

	# annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
#     annotation = read.delim(file.path(file.path("../example_data", "scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))	
    annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))
    
#     data_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", "Test_5_Example_scRNA_seq_Data")
    data_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", states_output_folder)
    
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
#     cell_types = c("Malignant")
#     cell_types = c("Bcell")

#     cell_types = c("Mast", "Monocyte", "Plasma", "Stellate")
#     cell_types = c("Mast", "Plasma", "Platelets", "Stellate")
    
	for(cell_type in cell_types)
	{
        
        # NEW CODE ADDED TO ACCOUNT FOR ONE STATE
#         input_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", input_folder_file_path, cell_type)
#         scaled_data = read.delim(file.path(input_dir, "expression_top_genes_scaled_filt.txt"))
#         annotationCT = read.delim(file.path(input_dir, "initial_state_assignment.txt"))
#         annotationCT = annotationCT[match(colnames(scaled_data), annotationCT$ID),]
        
#         annotationStates = annotationCT$State
#         annotationStatesNoDup = annotationStates[!duplicated(annotationStates)]
#         annotationStatesNoDupLen = length(annotationStatesNoDup)
#         ####
                                  
#         if (annotationStatesNoDupLen == 1) {
            
            
#             cat(paste0("\nConducting One State NMF for ", cell_type, "\n"))

#             scaledExpMat = read.delim(file.path(data_dir, cell_type, "expression_top_genes_scaled_step_5P1.txt"))
#             cellTypeBinaryH = read.delim(file.path(data_dir, cell_type, "Binary_H.txt"), row.names = 1)
#             avgGeneExpMatW = NMFGenerateW_One_State(cellTypeBinaryH, scaledExpMat)

#             write.table(avgGeneExpMatW, file.path(data_dir, cell_type, "Average_Cell_State_Gene_Expression_Matrix_W.txt"), sep = "\t")

        
        
#         } else {
        
            
#             cat(paste0("\nConducting Regular NMF for ", cell_type, "\n"))

#             scaledExpMat = read.delim(file.path(data_dir, cell_type, "expression_top_genes_scaled_step_5P1.txt"))
#             cellTypeBinaryH = read.delim(file.path(data_dir, cell_type, "Binary_H.txt"), row.names = 1)
#             avgGeneExpMatW = NMFGenerateW(cellTypeBinaryH, scaledExpMat)

#             write.table(avgGeneExpMatW, file.path(data_dir, cell_type, "Average_Cell_State_Gene_Expression_Matrix_W.txt"), sep = "\t")

        
#         }
        
        ### ADD NEW CODE FOR ERROR TESTING
#         previously_done_cell_types = c("Acinar", "Bcell", "Endothelial", "Fibroblast", "Malignant")
#         if (cell_type %in% previously_done_cell_types) {
#             cat(paste0("\nSkipping NMF for ", cell_type, "\n"))
#         } else {
#            cat(paste0("\nConducting Regular NMF for ", cell_type, "\n"))

#            scaledExpMat = read.delim(file.path(data_dir, cell_type, "expression_top_genes_scaled_step_5P1.txt"))
#            cellTypeBinaryH = read.delim(file.path(data_dir, cell_type, "Binary_H.txt"), row.names = 1)
#            avgGeneExpMatW = NMFGenerateW(cellTypeBinaryH, scaledExpMat)

#            write.table(avgGeneExpMatW, file.path(data_dir, cell_type, "Average_Cell_State_Gene_Expression_Matrix_W.txt"), sep = "\t")

#         }
        
        ################
        cat(paste0("\nConducting Regular NMF for ", cell_type, "\n"))

        # DOWNSIZING CONTRAINTS FOR LARGE MATRIX
#         scaledExpMat = read.delim(file.path(data_dir, cell_type, "expression_top_genes_scaled_step_5P1_DS.txt"), row.names = 1)
#         cellTypeBinaryH = read.delim(file.path(data_dir, cell_type, "Binary_H_DS.txt"), row.names = 1)

        scaledExpMat = read.delim(file.path(data_dir, cell_type, "expression_top_genes_scaled_step_5P1.txt"))
        cellTypeBinaryH = read.delim(file.path(data_dir, cell_type, "Binary_H.txt"), row.names = 1)
        avgGeneExpMatW = NMFGenerateW(cellTypeBinaryH, scaledExpMat)

        write.table(avgGeneExpMatW, file.path(data_dir, cell_type, "Average_Cell_State_Gene_Expression_Matrix_W.txt"), sep = "\t")


        # PushToJobQueue(paste("python lib/Python_Scripts/Test_Script.py", "arg1", "arg2", "arg3"))	
        
	}	
	RunJobQueue() 
	    
}else{
	cat("Skipping step 5 P5 (generating average gene expression matrix W)...\n")
}	

# if(!5 %in% skip_steps)
# {
# 	cat("\nStep 5 (cell state re-discovery in expression matrices)...\n")

# 	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# 	for(cell_type in key[,1])
# 	{	
# 		cat(paste("Extracting marker genes for cell states defined in:", cell_type, "\n"))
# 		n_clusters = key[key[,1] == cell_type, 2]
# 		PushToJobQueue(paste("Rscript state_discovery_extract_features_scRNA.R", discovery, fractions, cell_type, n_clusters)) 		 
# 	}	
# 	RunJobQueue()
	
# 	cat("\nStep 5 (cell state re-discovery in expression matrices): Running NMF on expression matrix...\n")

# 	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# 	for(cell_type in key[,1])
# 	{
# 		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, "expression_top_genes_scaled.txt")))
# 		{			
# 			next
# 		}
		
# 		n_clusters = key[key[,1] == cell_type, 2]
# 		for(restart in 1:nmf_restarts)
# 		{
# 			if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "restarts", restart, "estim.RData")))
# 			{
# 				PushToJobQueue(paste("Rscript state_discovery_NMF.R", "discovery", discovery, fractions, cell_type, n_clusters, restart))
# 			}else{					
# 				cat(paste0("Warning: Skipping NMF on '", cell_type, "' (number of states = ", n_clusters, ", restart ", restart, "), as the output file '", file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, n_clusters, "restarts", restart, "estim.RData"), "' already exists!\n"))
# 			}
# 		} 
# 	} 
# 	RunJobQueue()
		
# 	cat("Step 5 (cell state re-discovery in expression matrices): Aggregating NMF results...\n")
# 	for(cell_type in key[,1])
# 	{	
# 		if(!file.exists(file.path("../EcoTyper", discovery, fractions, "Cell_States", "discovery", cell_type, "expression_top_genes_scaled.txt")))
# 		{			
# 			next
# 		}				
# 		PushToJobQueue(paste("Rscript state_discovery_combine_NMF_restarts.R", "discovery", discovery, fractions, cell_type, max_clusters, nmf_restarts))
# 	} 
# 	RunJobQueue()
# 	cat("Step 5 (cell state re-discovery in expression matrices) finished successfully!\n")
# }else{
# 	cat("Skipping step 5 (cell state re-discovery in expression matrices)...\n")
# }	

# STEP 6

# DEFINE INPUT VARIABLES FOR GENERATE W FUNCTION

if(!6 %in% skip_steps) 
{
	# cat("\nStep 2 (cell state discovery on correrlation matrices): Calculating correlation matrices...\n")
    cat("\nStep 6 (generating initial gene info files)...\n")

	# annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
#     annotation = read.delim(file.path(file.path("../example_data", "scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))	
      annotation = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path))

#     data_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", "Test_5_Example_scRNA_seq_Data")
    
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
#     cell_types = c('Bcell', 'Malignant')
    
	for(cell_type in cell_types)
	{
        
        cat(paste0("\nGenerating initial gene info text files for ", cell_type, "\n"))
                
        PushToJobQueue(paste("Rscript state_discovery_initial_plots_Predefined_States_Mode.R", cell_type, states_output_folder, annotation_file_path))	
        
	}	
	RunJobQueue() 
	    
}else{
	cat("Skipping Step 6 (generating initial gene info files)...\n")
}	

# if(!6 %in% skip_steps) 
# {
# 	cat("\nStep 6 (extracting information for re-discovered cell states)...\n")

# 	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# 	for(cell_type in key[,1])
# 	{	
# 		cat(paste("Extracting cell states information for:", cell_type, "\n"))
# 		n_clusters = key[key[,1] == cell_type, 2]
# 		PushToJobQueue(paste("Rscript state_discovery_initial_plots.R", "discovery", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
# 	}	
# 	RunJobQueue()
# 	cat("Step 6 (extracting information for re-discovered cell states) finished successfully!\n")
# }else{
# 	cat("\nSkipping step 6 (extracting information for re-discovered cell states)...\n")
# }

# STEP 7

if(!7 %in% skip_steps)
{
	cat("\nStep 7 (cell state QC filter)...\n")

# 	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
        key = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", states_output_folder, "rank_data.txt"))
    
	for(cell_type in key[,1])
	{	
		cat(paste("Filtering low-quality cell states for:", cell_type, "\n"))
		n_clusters = key[key[,1] == cell_type, 2]
		PushToJobQueue(paste("Rscript state_discovery_first_filter_scRNA_Predefined_States_Mode.R", states_output_folder, annotation_file_path, cell_type, "State", paste(additional_columns, collapse = " "))) 		 
	}	
	RunJobQueue()
	cat("Step 7 (cell state QC filter) finished successfully!\n")
}else{
	cat("\nSkipping step 7 (cell state QC filter)...\n")
}

# STEP 7

# if(!7 %in% skip_steps)
# {
# 	cat("\nStep 7 (cell state QC filter)...\n")

# 	key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# 	for(cell_type in key[,1])
# 	{	
# 		cat(paste("Filtering low-quality cell states for:", cell_type, "\n"))
# 		n_clusters = key[key[,1] == cell_type, 2]
# 		PushToJobQueue(paste("Rscript state_discovery_first_filter_scRNA.R", discovery, fractions, cell_type, n_clusters, "State", paste(additional_columns, collapse = " "))) 		 
# 	}	
# 	RunJobQueue()
# 	cat("Step 7 (cell state QC filter) finished successfully!\n")
# }else{
# 	cat("\nSkipping step 7 (cell state QC filter)...\n")
# }


# STEP 8


if(!81 %in% skip_steps)
{
    	cat("\nStep 8 P1 (ecotype discovery)...\n")
# 	PushToJobQueue(paste("Rscript ecotypes_scRNA_Predefined_States_Mode.R", discovery, fractions, p_value_cutoff)) 
    
    #### OLD ECOTYPER PDSM
    
#     PushToJobQueue(paste("Rscript ecotypes_scRNA_Predefined_States_Mode.R", states_output_folder, annotation_file_path, p_val_cutoff, min_states)) 
    
    #### UPDATED ECOTYPER PDSM 
    
    PushToJobQueue(paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/ecotypes_scRNA_PDSM_Updated.R", states_output_folder, annotation_file_path, p_val_cutoff, min_states)) 
    
	RunJobQueue()
    cat("Step 8 P1 (ecotype discovery) finished successfully!\n")

}else{
	cat("\nSkipping step 8 P1 (ecotype discovery)...\n")
}

if(!82 %in% skip_steps)
{
    cat("Step 8 P2 (ecotype assignment) assigning ecotypes...!\n")
# 	PushToJobQueue(paste("Rscript ecotypes_assign_samples_scRNA_Predefined_States_Mode.R", discovery, fractions, "State",paste(additional_columns, collapse = " "))) 
    
    #### OLD ECOTYPER PDSM
    
#     PushToJobQueue(paste("Rscript ecotypes_assign_samples_scRNA_Predefined_States_Mode.R", states_output_folder, "State", annotation_file_path, paste(additional_columns, collapse = " "))) 
    
    #### UPDATED ECOTYPER PDSM
    
    PushToJobQueue(paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/ecotypes_assign_samples_scRNA_PDSM_Updated.R", states_output_folder, "State", annotation_file_path, paste(additional_columns, collapse = " "))) 
	
	RunJobQueue()
    cat("Step 8 P2 (ecotype assignment) ecotype assignment finished successfully!\n")

}else{
	cat("\nSkipping step 8 P2 (ecotype assignment)...\n")
}


# STEP 8

# if(!8 %in% skip_steps)
# {
# 	cat("\nStep 8 (ecotype discovery)...\n")
# # 	PushToJobQueue(paste("Rscript ecotypes_scRNA_Predefined_States_Mode.R", discovery, fractions, p_value_cutoff)) 
    
# #     PushToJobQueue(paste("Rscript ecotypes_scRNA_Predefined_States_Mode.R", states_output_folder, annotation_file_path, p_val_cutoff, min_states)) 
    
#     #### UPDATED ECOTYPER PDSM 
    
#     PushToJobQueue(paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/ecotypes_scRNA_PDSM_Updated.R", states_output_folder, annotation_file_path, p_val_cutoff, min_states)) 
    
# 	RunJobQueue()
#     cat("Step 8 (ecotype discovery) finished successfully!\n")
    
#     cat("Step 8 (ecotype discovery) assigning ecotypes...!\n")
# # 	PushToJobQueue(paste("Rscript ecotypes_assign_samples_scRNA_Predefined_States_Mode.R", discovery, fractions, "State",paste(additional_columns, collapse = " "))) 
    
# #     PushToJobQueue(paste("Rscript ecotypes_assign_samples_scRNA_Predefined_States_Mode.R", states_output_folder, "State", annotation_file_path, paste(additional_columns, collapse = " "))) 
    
#     #### UPDATED ECOTYPER PDSM
    
#     PushToJobQueue(paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/ecotypes_assign_samples_scRNA_PDSM_Updated.R", states_output_folder, "State", annotation_file_path, paste(additional_columns, collapse = " "))) 
	
# 	RunJobQueue()
#     cat("Step 8 (ecotype discovery) ecotype assignment finished successfully!\n")
    
# }else{
# 	cat("Skipping step 8 (ecotype discovery)...\n")
# }

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

end = Sys.time()
# cat(paste0("\nEcoTyper finished succesfully! Please find the results in: '", final_output, "'.\nRun time: ", format(end - start, digits = 0), "\n"))

