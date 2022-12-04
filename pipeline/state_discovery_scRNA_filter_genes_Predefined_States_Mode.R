suppressPackageStartupMessages({
library(data.table)
source("lib/misc.R")
source("lib/read_clinical.R")
})

# args = c("scRNA_CRC_Park", "scRNA_specific_genes", "B.cells") 
args = commandArgs(T)  
# dataset = args[1]
# fractions = args[2]
# cell_type = args[3]
cell_type = args[1]

annotation_file_path = args[2]
expression_matrix_file_path = args[3]
output_file_path = args[4]

# scaling_column = args[4] 
scaling_column = args[5]

# dataset_type = "discovery"
 
# input_dir = file.path("../datasets", dataset_type, dataset) 
# input_dir = file.path("../example_data") 
annotation_input_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", annotation_file_path)
raw_input_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", expression_matrix_file_path)


# output_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "Cell_type_specific_genes") 
# output_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", "Test_5_Example_scRNA_seq_Data", cell_type)
output_dir = file.path("../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", output_file_path, cell_type)

dir.create(output_dir, recursive = T, showWarning = F)

# annotation = read.delim(file.path(input_dir, paste0("annotation.txt")))
# annotation = read.delim(file.path(input_dir, paste0("scRNA_CRC_Annotation_Ecotype_States_Filtered.txt")))
annotation = read.delim(annotation_input_dir)

# raw_input = fread(file.path(input_dir, paste0("data.txt")), data.table  = F)
# raw_input = fread(file.path(input_dir, paste0("scRNA_CRC_data_Filtered.txt")), data.table  = F)
raw_input = fread(raw_input_dir, data.table  = F)

rownames(raw_input) =  raw_input[,1]
raw_input = raw_input[,-1]
annotation = annotation[annotation$ID %in% colnames(raw_input),]

tb = as.data.frame(table(as.character(annotation$CellType)))

# for(i in 1:nrow(tb))
# {
# 	if(tb[i,2] < 50)
# 	{
# 		warning(paste0("Only ", tb[i,2], " single cells are available for cell type: ", tb[i,1], ". At least 50 are required. Skipping this cell type from the EcoTyper analysis!"))
# 	}
# }

# tb = tb[tb[,2] >= 50,]
# if(nrow(tb) == 0)
# {
# 	stop("No cell type has at least 50 cells! EcoTyper cannot run!")
# }

annotation = annotation[annotation$CellType %in% tb[,1],]
splits = split(annotation, annotation$CellType)

set.seed(1234)
annotation = do.call(rbind, lapply(splits, function(spl){
	spl[sample(1:nrow(spl), min(500, nrow(spl))),]
	}))

raw_input = raw_input[,colnames(raw_input) %in% annotation$ID]
#raw_input = raw_input[apply(raw_input,  1, var) > 0 && apply(raw_input, 1, function(x) sum(x > 0) >= 3),] 

log_data = log2(raw_input + 1)

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
gene_info = doDE(scaled_data, clinical$CellType, cell_type)

colnames(gene_info)[2] = "CellType"
write.table(gene_info, file.path(output_dir, paste0(cell_type, "_cell_type_specific_genes_raw.txt")), sep = "\t", row.names = F)
gene_info = gene_info[gene_info$Q <= 0.05 & gene_info$FC > 0,]
write.table(gene_info, file.path(output_dir, paste0(cell_type, "_cell_type_specific_genes.txt")), sep = "\t", row.names = F)
