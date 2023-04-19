suppressPackageStartupMessages({
library(argparse)
source("pipeline/lib/config.R") 
source("pipeline/lib/multithreading.R")
})

parser <- ArgumentParser(add_help = F)

arguments = parser$add_argument_group('Arguments')

# arguments$add_argument("-d", "--discovery", type="character", default="Carcinoma", 
#     help=paste0("The name of the discovery dataset used to define cell states and ecotypes. Accepted values: ",
#     "'Carcinoma' will recover the cell states and ecotypes defined across carcinomas, as described in the EcoTyper carcinoma paper, ",
#     "'Lymphoma' will recover the cell states and ecotypes defined in diffuse large B cell lymphoma (DLBCL), as described in the EcoTyper lymphoma paper, ",
#     "'<MyDiscovery>' the value used in the field 'Discovery dataset name' of the config file used for running EcoTyper discovery ('EcoTyper_discovery.R') script. ",
#     "[default: '%(default)s']"),
#     metavar="<character>") 
arguments$add_argument("-m", "--matrix", type = "character", metavar="<PATH>",
    help="Path to a tab-delimited file containing the input bulk tissue expression matrix, with gene names on the first column and sample ids as column names [required].")
arguments$add_argument("-a", "--annotation", type = "character",  metavar="<PATH>", default="NULL",
    help="Path to a tab-delimited annotation file containing the annotation of samples in the input matrix. This file has to contain in column 'ID' the same ids used as column names in the input matrix, and any number of additional columns. The additional columns can be plotted as color bars in the output heatmaps. [default: '%(default)s']")
# arguments$add_argument("-c", "--columns",  type = "character",  metavar="<character>", default="NULL", 
#     help="A comma-spearated list of column names from the annotation file to be plotted as color bar in the output heatmaps. [default: '%(default)s']")
arguments$add_argument("-t", "--threads",  type = "integer",  metavar="<integer>", default=10, 
    help="Number of threads. [default: '%(default)s']")
# arguments$add_argument("-o", "--output", type = "character",  metavar="<PATH>", default="RecoveryOutput",
#     help="Output directory path. [default: '%(default)s']")
# arguments$add_argument("-h", "--help", action='store_true', help="Print help message.")


#### ADD NEW ARGUMENTS
arguments$add_argument("-bi", "--bulkInputsDir", type = "character", metavar="<PATH>")
arguments$add_argument("-dm", "--discoveryMatrix", type = "character", metavar="<PATH>")
arguments$add_argument("-da", "--discoveryAnnotation", type = "character", metavar="<PATH>")
arguments$add_argument("-k", "--keyDir", type = "character", metavar="<PATH>")
arguments$add_argument("-o", "--outputDir", type = "character",  metavar="<PATH>", default="RecoveryOutput")
arguments$add_argument("-s", "--skipStep", type = "integer",  metavar="<integer>", default=0)
arguments$add_argument("-pd", "--priorDisc", type = "character",  metavar="<PATH>", default="DiscoveryOutput")
arguments$add_argument("-pr", "--priorRank", type = "character",  metavar="<PATH>", default="DiscoveryOutput")
arguments$add_argument("-ssm", "--singleStateMode", type = "integer",  metavar="<integer>", default=0)


args <- parser$parse_args()
#print(args)
if(is.null(args$matrix))
{
	parser$print_help()
	quit()
}

bulkInputsDir = args$bulkInputsDir
discoveryMatrix = args$discoveryMatrix
discoveryAnnotation = args$discoveryAnnotation
keyDir = args$keyDir
outputDir = args$outputDir
skipStep = args$skipStep
priorDisc = args$priorDisc
priorRank = args$priorRank
singleStateMode = args$singleStateMode

# annPath = "Generate_Cell_State_Abundances_Predefined_States_Test_Outputs/Test_8_Peng_Discovered_States_EcoTyper_PDSM_Run_083022/Single_Cell_Inputs/scRNA_Seq_Peng_Discovered_States_Annotations.txt"

# inputMatPath = "Generate_Cell_State_Abundances_Predefined_States_Test_Outputs/Test_8_Peng_Discovered_States_EcoTyper_PDSM_Run_083022/Single_Cell_Inputs/scRNA_Seq_Peng_Discovered_States_Count_Matrix.txt"

input_mat = normalizePath(args$matrix)
# input_mat = normalizePath("tcga_PDSM.txt")

# input_mat = normalizePath(inputMatPath)

# discovery = args$discovery

annotation_path = normalizePath(args$annotation)
# annotation_path = normalizePath("tcga_annotations_PDSM.txt")
# annotation_path = normalizePath(annPath)

# columns = args$columns

n_threads = args$threads
# n_threads = 8

# if(!file.exists(input_mat))
# {
# 	stop("Error: Path to the input expression matrix does not exist!")
# }

# if(! discovery %in% c("Carcinoma", "Lymphoma"))
# {
# 	config_file = file.path("EcoTyper", discovery, "config_used.yml")	
# 	if(!file.exists(config_file))
# 	{
# 		stop("Error: Cannot read the config file used for the discovery of cell states and ecotypes. It should be in the following path: '", config_file, "'. Please make sure that the '--discovery (-d)' argument provided is correct!")
# 	}
# 	config <- config::get(file = config_file)
	
# 	fractions = config$"Input"$"Cell type fractions" 
# 	if(is.null(fractions))
# 	{
# 		if(config$"Pipeline settings"$"Filter non cell type specific genes")
# 		{
# 			fractions = "Cell_type_specific_genes"
# 		}else{
# 			fractions = "All_genes"
# 		}
# 	}else{
# 		if(!fractions %in% c("Carcinoma_Fractions", "Lymphoma_Fractions"))
# 		{
# 			fractions = "Custom"
# 		}
# 	}
# }else{
# 	fractions = paste0(discovery, "_Fractions") 
# }

# if(annotation_path != "NULL" && columns != "NULL")
# {
# 	if(!file.exists(annotation_path))
# 	{
# 		stop("Error: Path to the input annotation file does not exist!")
# 	}else{
# 		annotation = read.delim(annotation_path)
# 		#print(head(annotation))
# 		if(!"ID" %in% colnames(annotation))
# 		{
# 			stop("Error: The annotation file provided, does not contain the column 'ID'.")
# 		}

# 		additional_columns = strsplit(columns, ",")[[1]]
		
# 		if(!all(additional_columns %in% colnames(annotation)))
# 		{
# 			stop(paste0("The following columns are missing from the annotation file provided: ", "'", 
# 				paste(additional_columns[!additional_columns %in% colnames(annotation)], collapse = "'"), "'"))			
# 		} 
# 	}
# }else{
# 	additional_columns = c()
# }

# recovery = gsub(".tsv$", "", gsub(".txt$", "", basename(input_mat)))

# input_dir = file.path("datasets", "bulk", recovery)
# input_dir = file.path("Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", "Test_8_Peng_Discovered_States_EcoTyper_PDSM_Run_083022", "TCGA_Recovery_090122", "Inputs")
input_dir = file.path("Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", args$bulkInputsDir)

# dir.create(input_dir, recursive = T, showWarning = F)

# dir.create(file.path(args$output, recovery), recursive = T, showWarnings = F)
# final_output = normalizePath(file.path(args$output, recovery))

# PushToJobQueue(paste0("ln -sf ", input_mat, " ", file.path(input_dir, "data.txt")))
PushToJobQueue(paste0("ln -sf ", input_mat, " ", file.path(input_dir, "data.txt")))
RunJobQueue()


# PushToJobQueue(paste0("ln -sf ", annotation_path, " ", file.path(input_dir, "annotation.txt")))
PushToJobQueue(paste0("ln -sf ", annotation_path, " ", file.path(input_dir, "annotation.txt")))
RunJobQueue()

# start = Sys.time()
cur_dir = getwd()
setwd("pipeline")

# key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# key = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", "Test_8_Peng_Discovered_States_EcoTyper_PDSM_Run_083022", "rank_data.txt"))

key = read.delim(file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", keyDir, "rank_data.txt"))

#### single state cell types created manually - automize this later on

if (singleStateMode == 0) {
    singleStateCellTypes =c("No_SSM")
    
} else {
    singleStateCellTypes = c("Mast", "Plasma", "Stellate", "Endothelial", "Acinar", "Normal_Epithelial", "B_cell")   

}
    

priorKey = read.delim(file.path(priorRank, "rank_data.txt"))

if(skipStep != 1) {
    
    for(cell_type in key[,1])
    {
        #### NEW CODE ADDED TO ACCOUNT FOR ONE STATE
        
        if(cell_type %in% singleStateCellTypes) {
            
            n_states = priorKey[priorKey[,1] == cell_type, 2]
            PushToJobQueue((paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/state_recovery_bulk_PDSM_Updated_Single_State.R", cell_type, keyDir, bulkInputsDir, outputDir, priorDisc, n_states))) 
            
        } else {
            
          n_states = key[key[,1] == cell_type, 2]
        #     PushToJobQueue((paste("Rscript state_recovery_bulk.R", discovery, fractions, cell_type, n_states, recovery, "FALSE", paste(additional_columns, collapse = " ")))) 
        
        #### OLD STATE PDSM RECOVERY 
#         PushToJobQueue((paste("Rscript state_recovery_bulk_Predefined_States_Mode.R", cell_type, keyDir, bulkInputsDir, outputDir))) 
    
        #### UPDATED STATE PDSM RECOVERY
        PushToJobQueue((paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/state_recovery_bulk_PDSM_Updated.R", cell_type, keyDir, bulkInputsDir, outputDir))) 

        
        }
        
    }   
    RunJobQueue()
    cat("\nFinished state recovery on bulk dataset, moving onto ecotype recovery\n")
}else{
	cat("Skipping state recovery on bulk dataset\n")
}	


# PushToJobQueue((paste("Rscript ecotypes_recovery.R", discovery, fractions, recovery, paste(additional_columns, collapse = " ")))) 
# RunJobQueue()

if(skipStep != 2) {

          #### OLD ECOTYPE PDSM RECOVERY 
#           PushToJobQueue((paste("Rscript ecotypes_recovery_Predefined_States_Mode.R", keyDir, bulkInputsDir, outputDir))) 
#           RunJobQueue()
#         #### UPDATED ECOTYPE PDSM RECOVERY
        PushToJobQueue((paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/ecotypes_recovery_PDSM_Updated.R", keyDir, bulkInputsDir, outputDir))) 
        RunJobQueue()
    
    cat("\nEcotype recovery finished successfully\n")
#     if (singleStateMode == 0) {
        
            
#         cat("\nPerforming ecotype recovery on bulk dataset\n")
#         #### OLD ECOTYPE PDSM RECOVERY 
#     #     PushToJobQueue((paste("Rscript ecotypes_recovery_Predefined_States_Mode.R", keyDir, bulkInputsDir, outputDir))) 

#         #### UPDATED ECOTYPE PDSM RECOVERY
#         PushToJobQueue((paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/ecotypes_recovery_PDSM_Updated.R", keyDir, bulkInputsDir, outputDir))) 

#         RunJobQueue()

#         # cat("\nCopying EcoTyper results to the output folder!\n")
#         cat("\nEcotype recovery finished successfully\n")
    
#     } else {

#         cat("\nConducting ecotype recovery with single states\n")
        
#         PushToJobQueue((paste("Rscript /duo4/users/pchati/EcoTyper_Updated/ecotyper/pipeline_PDSM/ecotypes_recovery_PDSM_Updated_Single_State.R", keyDir, bulkInputsDir, outputDir, priorDisc, priorRank, singleStateCellTypes))) 
        
#         RunJobQueue()
        
#         cat("\nEcotype recovery single states mode finished successfully\n")
#     }


}else{
	cat("Skipping ecotype recovery on bulk dataset\n")
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

# key = read.delim(file.path("../EcoTyper", discovery, fractions, "Analysis", "rank_selection", "rank_data.txt"))
# for(cell_type in key[,1])
# {		
# 	n_clusters = key[key[,1] == cell_type, 2]	
# 	ct_output = file.path(final_output, cell_type)
# 	dir.create(ct_output, recursive = T, showWarnings = F)
	
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "state_abundances.txt"), ct_output))
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "state_assignment.txt"), ct_output))
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "state_assignment_heatmap.pdf"), ct_output))	
# 	system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Cell_States", "recovery", recovery, cell_type, n_clusters, "state_assignment_heatmap.png"), ct_output))	
# }	

# ct_output = file.path(final_output, "Ecotypes")
# dir.create(ct_output, recursive = T, showWarnings = F)
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "ecotype_assignment.txt"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "ecotype_abundance.txt"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "heatmap_assigned_samples_viridis.pdf"), ct_output))
# system(paste("cp -f", file.path("../EcoTyper", discovery, fractions, "Ecotypes", "recovery", recovery, "heatmap_assigned_samples_viridis.png"), ct_output))

# end = Sys.time()
# cat(paste0("\nEcoTyper finished succesfully! Please find the results in: '", final_output, "'.\nRun time: ", format(end - start, digits = 0), "\n"))

