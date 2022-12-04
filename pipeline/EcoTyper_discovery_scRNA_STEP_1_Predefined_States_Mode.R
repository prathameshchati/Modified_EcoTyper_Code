suppressPackageStartupMessages({
library(config)
library(argparse)
source("pipeline/lib/config.R") 
source("pipeline/lib/misc.R") 
source("pipeline/lib/multithreading.R")
source("pipeline/lib/ecotyper_NMF_Generate_W_Function.R")
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
check_discovery_configuration_scRNA(config)

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

if(!1 %in% skip_steps & config$"Pipeline settings"$"Filter non cell type specific genes")
{
	cat("\nStep 1 (extract cell type specific genes)...\n")

	# annotation = read.delim(file.path(file.path("../datasets/discovery", discovery, "annotation.txt")))	
    annotation = read.delim(file.path(file.path("../example_data", "scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))	
	cell_types = unlist(levels(as.factor(as.character(annotation$CellType))))	
	for(cell_type in cell_types)
	{
		print(cell_type)
		# PushToJobQueue(paste("Rscript state_discovery_scRNA_filter_genes.R", discovery, fractions, cell_type, scale_column))	
        PushToJobQueue(paste("Rscript state_discovery_scRNA_filter_genes_Predefined_States_Mode.R", discovery, fractions, cell_type, scale_column))	
        
	}
	RunJobQueue()	
	cat("Step 1 (extract cell type specific genes) finished successfully!\n")
	
}else{
	cat("Skipping step 1 (extract cell type specific genes)...\n")
}
