suppressPackageStartupMessages({
library(cluster)
library(ggplot2)
library(viridis)
library(reshape2)
source("lib/misc.R")
source("lib/heatmaps.R")
source("lib/read_clinical.R")
})

# args = c("discovery_scRNA_CRC", "Cell_type_specific_genes", "Ecotype", "Tissue")
args = commandArgs(T) 
# dataset = args[1]
# fractions = args[2]
# top_cols = args[3:length(args)]

states_output_folder = args[1]
fraction_processing = args[2]
annotation_file_path = args[3]
top_cols = args[4:length(args)]

if(is.na(top_cols[1]))
{
	top_cols = c("Ecotype")
}
top_cols = unique(c(top_cols, "Ecotype"))

# key_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "rank_selection")
# states_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery")
# output_dir = file.path("../EcoTyper", dataset, fractions, "Ecotypes", "discovery")

key_dir = file.path("/duo4/users/pchati/ecotyper/Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", states_output_folder)
states_dir = file.path("/duo4/users/pchati/ecotyper/Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", states_output_folder)
output_dir = file.path("/duo4/users/pchati/ecotyper/Generate_Cell_State_Abundances_Predefined_States_Test_Outputs", states_output_folder, "Ecotypes")

dir.create(output_dir, recursive = T, showWarning = F) 

key = read.delim(file.path(key_dir, "rank_data.txt"))
ecotypes = read.delim(file.path(output_dir, 'ecotypes.txt'), sep = '\t')
mapping = read.delim(file.path(output_dir, 'mapping_all_states.txt'), sep = '\t')
ecotypes$Ecotype = ecotype_to_factor(ecotypes$Ecotype)
ecotypes = ecotypes[order(ecotypes$Ecotype),]

all_H = NULL
all_classes_filt = NULL
cell_ids = c()
for(cell_type in key[,1])
{	
	#print(cell_type) 
	n_clusters = key[key[,1] == cell_type,2]
    
# 	mapping = read.delim(file.path(states_dir, cell_type, n_clusters, 'mapping_to_initial_states.txt'))
# 	classes = read.delim(file.path(states_dir, cell_type, n_clusters, 'initial_state_assignment.txt'))
    
	mapping = read.delim(file.path(states_dir, cell_type, 'mapping_to_initial_states.txt'))
	classes = read.delim(file.path(states_dir, cell_type, 'initial_state_assignment.txt'))

# 	clinical = read_clinical(classes$ID, dataset = dataset)
    clinical = read_clinical_manual(classes$ID, annotation_file_path = annotation_file_path, dataset = dataset)
    
	classes$Sample = clinical$Sample
	
	cell_ids = unique(c(cell_ids, classes$ID))

	classes = as.data.frame(table(classes$Sample, classes$State))
	colnames(classes) = c("ID", "State", "Freq")
	splits = split(classes, classes$ID)
	classes = do.call(rbind, lapply(splits, function(spl)
	{
		spl$Frac = spl$Freq / sum(spl$Freq)
		spl$Max= ifelse(which.max(spl$Freq) == spl$Freq, 1, 0)
		spl
		}))
	classes$CellType = cell_type		
	H = dcast(classes, State~ID, value.var = "Frac")
	rownames(H)  = H[,1]
	H = H[,-1]
	H_raw = H
	H =H[match(mapping$InitialState, rownames(H)),]
	rownames(H) = mapping$State

	rownames(H) = paste0(cell_type, "_", rownames(H))	
	keep_rowname = c(rownames(all_H), rownames(H))
	all_H = rbind.fill(all_H, H)
	rownames(all_H) = keep_rowname

	classes = as.data.frame(apply(H_raw, 2, function(x) {
		idx = which.max(x)
		if(length(idx)== 0)
		{
			"Unassigned"
		}else{
			rownames(H_raw)[idx]
		}
		
	}
		))
	classes_raw = data.frame(ID = rownames(classes), InitialState = classes[,1])
	classes_raw = classes_raw[classes_raw$InitialState %in% mapping$InitialState,]
	classes_raw$State = mapping[match(classes_raw$InitialState, mapping$InitialState), "State"]
	classes = classes_raw
	classes = classes[,c("ID", "State")]
	clusters = ecotypes[ecotypes$CellType == cell_type,] 
	classes = classes[classes$State %in% clusters$State,]
	
	colnames(classes) = c('ID', cell_type)

	if(is.null(all_classes_filt))
	{
		all_classes_filt = classes
	}else{
		all_classes_filt = merge(all_classes_filt, classes, by = 'ID', all = T)
	}
} 

all_H = all_H[match(ecotypes$ID, rownames(all_H)),]
write.table(all_H, file.path(output_dir, "combined_state_abundances.txt"), sep = "\t")

H = do.call(rbind, lapply(levels(ecotypes$Ecotype), function(clst){
	clst <<- clst
	#print(clst)
	inc <<- ecotypes[ecotypes$Ecotype == clst,]$ID
	apply(all_H[rownames(all_H) %in% inc,,drop = F], 2, mean, na.rm = T)
}))

rownames(H) = levels(ecotypes$Ecotype)
#write.table(H, file.path(output_dir, "ecotype_abundance.txt"), sep = "\t")
H = apply(H, 2, function(x) x / sum(x, na.rm = T))
write.table(H, file.path(output_dir, "ecotype_abundance.txt"), sep = "\t")

p_vals = do.call(rbind, lapply(levels(ecotypes$Ecotype), function(clst){
	clst <<- clst
	#print(clst)
	inc <<- ecotypes[ecotypes$Ecotype == clst,]$ID
	
	apply(all_H, 2, function(x){
		x <<- x
		err <<- F
		p <<- 1
		tryCatch({
			p <<- t.test(x[rownames(all_H) %in% inc], x[!(rownames(all_H) %in% inc)])$p.value
		}, error = function(x) err <<- T)
		if(err)
		{
			p = NA 
		}
		p
	})
	
})) 
rownames(p_vals) = levels(ecotypes$Ecotype)
write.table(p_vals, file.path(output_dir, "assignment_p_vals.txt"), sep = "\t")

assignment = as.data.frame(apply(H, 2, function(x) {
	idx = which.max(x)
	if(length(idx)==0)
	{
		"Unassigned"
	}else{
		rownames(H)[idx]
	}	
}))

clinical = data.frame(ID = rownames(assignment), MaxEcotype = assignment[,1])
clinical$AssignmentP = sapply(1:ncol(H), function(i) {
	idx = which.max(H[,i])
	if(length(idx)==0)
	{
		NA
	}else{
		p_vals[idx, i] 
	}		
})
clinical$AssignmentQ = p.adjust(clinical$AssignmentP, method = "BH")

clinical$AssignedToEcotypeStates = clinical$ID %in% all_classes_filt$ID
#clinical$Ecotype = ifelse((clinical$AssignmentQ < 0.25) & clinical$AssignedToEcotypeStates, as.character(clinical$MaxEcotype), "Unassigned")
clinical$Ecotype = ifelse( clinical$AssignedToEcotypeStates, as.character(clinical$MaxEcotype), "Unassigned")
clinical$Ecotype = factor(as.character(clinical$Ecotype), levels = c(levels(ecotypes$Ecotype), "Unassigned"))

# annotation = read_clinical(cell_ids, dataset = dataset)
annotation = read_clinical_manual(cell_ids, annotation_file_path = annotation_file_path, dataset = dataset)     
    
additional_columns = top_cols[top_cols %in% colnames(annotation)]
agg_ann = do.call(cbind, sapply(additional_columns, function(x){
	tb = as.data.frame(table(annotation$Sample, annotation[,x]))
	tb = tb[tb[,3] > 0,]
	if(sum(duplicated(tb[,1])) > 0){
		warning(paste0("Not plotting column '", x, "', as it has multiple distinct values within the same sample!"))
		return(NULL)
	}
	tb = tb[match(clinical$ID, tb[,1]),]
	as.character(tb[,2])
	}, simplify = F))
if(!is.null(agg_ann))
{
	clinical = cbind(clinical, agg_ann)
	top_cols = unique(c(colnames(agg_ann), "Ecotype"))	
}else{
	top_cols = "Ecotype"
}

#tmp = read_clinical(clinical$ID, dataset = dataset)
#to_rem = colnames(tmp)[colnames(tmp) %in% colnames(clinical)]
#to_rem = to_rem[to_rem != "ID"]
#tmp = tmp[,!colnames(tmp) %in% to_rem]
#clinical = merge(clinical, tmp, by = "ID", all.x = T)

clinical = clinical[order(clinical$Ecotype),]
H = H[,match(clinical$ID, colnames(H))]
all_H = all_H[,match(clinical$ID, colnames(all_H))]
all_H = all_H[match(ecotypes$ID, rownames(all_H)),]

write.table(clinical, file.path(output_dir, "initial_ecotype_assignment.txt"), sep = "\t")

rownames(clinical) = clinical$ID
rownames(ecotypes) = ecotypes$ID
    
#### NEW CODE TO AVOID NA VALUES IN MATRIX
all_H[is.na(all_H)] = 0

h <- heatmap_simple(all_H, top_annotation = clinical, top_columns = top_cols, 
	left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State"),
	column_split = ifelse(clinical$Ecotype == "Unassigned", "Unassigned", "Assigned"),
	width = unit(7, "in"), height = unit(4, "in"),
	legend_name = "State abundance",
	color_range = seq(0, quantile(as.matrix(all_H), .9, na.rm = T), length.out = 8), color_palette = c("gray", viridis(8)), raster_quality = 5)

pdf(file.path(output_dir, "heatmap_all_samples.pdf"), width = 12, height = 9)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
tmp = dev.off()

clinical_filt = clinical[clinical$Ecotype != "Unassigned",]
clinical_filt$Ecotype = factor(as.character(clinical_filt$Ecotype), levels = levels(ecotypes$Ecotype))
clinical_filt = clinical_filt[order(clinical_filt$Ecotype),]
write.table(clinical_filt, file.path(output_dir, "ecotype_assignment.txt"), sep = "\t")
small_H = as.matrix(all_H[,match(clinical_filt$ID, colnames(all_H))])
    
    
#### NEW CODE TO AVOID NA VALUES IN MATRIX
small_H[is.na(small_H)] = 0    

if(is.null(small_H) || nrow(small_H) == 0 || ncol(small_H) == 0)
{
	stop(paste("No samples were assigned to ecotypes!"))
}
h = heatmap_simple(small_H, top_annotation = clinical_filt, top_columns = top_cols, 
	left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State"),
	width = unit(5, "in"), height = unit(3, "in"),
	legend_name = "State abundance",
	color_range = seq(0, quantile(as.matrix(all_H), .9, na.rm = T), length.out = 8), color_palette = c("gray", viridis(8)), raster_quality = 20)

pdf(file.path(output_dir, "heatmap_assigned_samples_viridis.pdf"), width = 8, height = 8)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
#draw(h1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
suppressWarnings({
rect = rectangle_annotation_coordinates(ecotypes$Ecotype, clinical_filt$Ecotype)
})
decorate_heatmap_body("hmap", {
    grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3))
})
tmp = dev.off()

png(file.path(output_dir, "heatmap_assigned_samples_viridis.png"), width = 8, height = 8, units = "in", res = 300)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
	grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3)) 
})
tmp = dev.off()

h = heatmap_simple(small_H, top_annotation = clinical_filt, top_columns = top_cols, 
	left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State"),
	width = unit(5, "in"), height = unit(3, "in"), 
	legend_name = "State abundance",
	color_range = c(seq(0, quantile(small_H, .8, na.rm = T), length.out = 8)), color_palette = c("gray", brewer.pal(8, "YlGnBu")), raster_quality = 20)

pdf(file.path(output_dir, "heatmap_assigned_samples_YlGnBu.pdf"), width = 8, height = 6)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
    grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3))
})
tmp = dev.off()

png(file.path(output_dir, "heatmap_assigned_samples_YlGnBu.png"), width = 8, height = 6, units = "in", res = 300)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
    grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3))
})
tmp = dev.off()


