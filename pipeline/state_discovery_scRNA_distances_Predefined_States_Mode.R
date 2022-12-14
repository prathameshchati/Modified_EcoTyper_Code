suppressPackageStartupMessages({
library(data.table)
library(HiClimR)
source("lib/misc.R")
source("lib/read_clinical.R")
})

# args = c("scRNA_CRC_Park", "scRNA_specific_genes", "B.cells", "TRUE") 
args = commandArgs(T)  
# dataset = args[1]
# fractions = args[2]

# cell_type = args[3]
cell_type = args[1]

# filter_genes = as.logical(args[4])
filter_genes = as.logical(args[2])

annotation_file_path = args[3]
expression_matrix_file_path = args[4]
input_folder_file_path = args[5]
output_folder_file_path = args[6]

# scaling_column = args[5] 
scaling_column = args[7]

dataset_type = "discovery"
 
# dataset_dir = file.path("../datasets/discovery", dataset) 
# dataset_dir = file.path("../example_data") 
annotation_input_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path)
raw_input_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", expression_matrix_file_path)

# input_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "Cell_type_specific_genes")
# input_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", "Test_5_Example_scRNA_seq_Data", cell_type)
input_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", input_folder_file_path, cell_type)

# output_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery_cross_cor", cell_type)
# output_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", "Test_5_Example_scRNA_seq_Data", cell_type)
output_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", output_folder_file_path, cell_type)

# dir.create(output_dir, recursive = T, showWarning = F)

# annotation = read.delim(file.path(dataset_dir, paste0("annotation.txt")))
# annotation = read.delim(file.path(dataset_dir, paste0("scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))
annotation = read.delim(annotation_input_dir)

annotation = annotation[annotation$CellType == cell_type,]

# if(nrow(annotation) < 50)
# {
# 	stop(paste0("Only ", nrow(annotation), " single cells are available for cell type: ", cell_type, ". At least 50 are required. Skipping this cell type!\n"))	
# }

# if(nrow(annotation) > 2500)
# {
# 	warning(paste0("There are more than 2500 single cells available for cell type '", cell_type, "'. ", "Subsampling to 2500 cells.\n"))
# 	set.seed(1234)
# 	annotation = annotation[sample(1:nrow(annotation), 2500),]
# }

# raw_input = fread(file.path(dataset_dir, paste0("data.txt")), data.table  = F)
# raw_input = fread(file.path(dataset_dir, paste0("scRNA_CRC_data_Filtered.txt")), data.table  = F)
raw_input = fread(raw_input_dir, data.table  = F)

rownames(raw_input) =  raw_input[,1]
raw_input = raw_input[,-1]
raw_input = raw_input[,colnames(raw_input) %in% annotation$ID]

log_data = log2(raw_input + 1)

write.table(log_data, file.path(output_dir, "expression_full_matrix_log2.txt"), sep = "\t", row.names = T)

# clinical = read_clinical(colnames(log_data), dataset = dataset, dataset_type = "discovery")
# clinical = read_clinical_manual(colnames(log_data), dataset = dataset, dataset_type = "discovery")
clinical = read_clinical_manual(colnames(log_data), annotation_file_path = annotation_file_path, dataset = dataset, dataset_type = "discovery")

if(is.na(scaling_column)) 
{
	by = NULL
}else{
	by = as.character(clinical[,scaling_column])
}

scaled_data = scale_data(log_data, by = by)
scaled_data[is.na(scaled_data)] = 0
write.table(scaled_data, file.path(output_dir, "expression_full_matrix_scaled.txt"), sep = "\t")

if(filter_genes)
{
	cat(paste0("Filtering '", cell_type, "' profiles for cell type specific genes...\n")) 
	if(!file.exists(file.path(input_dir, paste0(cell_type, "_cell_type_specific_genes.txt"))))
	{
		stop(paste0("List of cell type specific genes not found for cell type '", cell_type, "'. Please make sure that step 1 finished sucessfully!\n"))
	}
	gene_info = read.delim(file.path(input_dir, paste0(cell_type, "_cell_type_specific_genes.txt")))
	gene_info = gene_info[gene_info$CellType == cell_type,]
	scaled_data = scaled_data[rownames(scaled_data) %in% gene_info$Gene,]
}else{
	cat(paste0("Not filtering '", cell_type, "' profiles for cell type specific genes, using the full transcriptome...\n"))
}

write.table(scaled_data, file.path(output_dir, "expression_top_genes_scaled_filt.txt"), sep = "\t")

# d = fastCor(scaled_data)
# write.table(d, file.path(output_dir, "expression_top_genes_scaled.txt"), sep = "\t")

