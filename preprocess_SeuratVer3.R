library(monocle)
library(FNN)
library(plyr)
library(Matrix)
library(Seurat)
source('/Users/fanny/Documents/Shendure/mouse/helper_functions.R')
setwd('/Users/fanny/Documents/Shendure/mouse')
setwd('/Users/fanny/Documents/Shendure/mouse/tabulamuris')
#setwd('/Users/fanny/Documents/Shendure/mouse/park')


# Load in the ATAC data (once)
activity_score_matrix = '/Users/fanny/Documents/Shendure/mouse/aggregate_counts.all_reads.rds'
activity_score_matrix = '/Users/fanny/Documents/Shendure/mouse/activity_scores.quantitative.rds'
atac.metadata.file = "/Users/fanny/Documents/Shendure/mouse/final_cell_type_assignments.txt"
atac.gene_activity_scores = readRDS(activity_score_matrix)
atac.metadata = as.data.frame(readr::read_delim(atac.metadata.file, '\t'))
rownames(atac.metadata) = atac.metadata$cell

# Load in the RNA data
rna.cds_list = readRDS('/Users/fanny/Documents/Shendure/rna.cds_list.microwell_seq.rds')
#rna.cds_list = readRDS('/Users/fanny/Documents/Shendure/mouse/tabulamuris/rna.cds_list.rds')
#rna.cds_list = readRDS('/Users/fanny/Documents/Shendure/mouse/park/rna.kidney.cds_list.rds')
tissues = names(rna.cds_list)

#filter low depth cells
atac = readRDS('/Users/fanny/Documents/Shendure/mouse/atac_matrix.binary.qc_filtered.rds')
atac.totals = Matrix::colSums(atac)

atac.depth_filter = 1800 #TODO:2000, 1800
rna.depth_filter = 600 #TODO:800, 600

atac_cells = colnames(atac)[atac.totals > atac.depth_filter]

rna_cells = ""
for (t in tissues){
  rna.temp = rna.cds_list[[t]]
  rna.totals.temp = Matrix::colSums(exprs(rna.temp))
  rna.temp.high_depth_cells = colnames(rna.temp)[rna.totals.temp > rna.depth_filter]
  rna_cells = c(rna_cells, rna.temp.high_depth_cells)
}

knn_k = 20 #TODO
cell_type_filter = 0.5 #TODO
tissues
#"Kidney", "Marrow", "Spleen", "Thymus", "Liver", "Heart", "Lung", "Brain", "Large_Intestine"
tissue = tissues[2]
# Get RNA data for selected tissue
rna.cds = rna.cds_list[[tissue]]
rna.cds = rna.cds[, !is.na(rna.cds$cell_type)]
rna.cds = rna.cds[,(rna.cds$cell %in% rna_cells)]

print(paste0('Processing: ', tissue, '...'))
dir.create(paste0('mlproject/'), showWarnings = TRUE)
dir.create(paste0('mlproject/', tissue), showWarnings = TRUE)

# Get ATAC data for selected tissue
if (tissue == 'Testis'){
  subset_meta = atac.metadata[grepl('Testes', atac.metadata$tissue),]
} else if (tissue == 'Marrow') {
  subset_meta = atac.metadata[grepl('BoneMarrow', atac.metadata$tissue),]
} else if (tissue == 'Large_Intestine') {
  subset_meta = atac.metadata[grepl('Intestine', atac.metadata$tissue),]
} else {
  subset_meta = atac.metadata[grepl(tissue, atac.metadata$tissue),]
}
atac.gene_activity_scores.sub = atac.gene_activity_scores[, subset_meta$cell]
atac.pd = new("AnnotatedDataFrame", data = subset_meta)
atac.fd = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(atac.gene_activity_scores.sub), gene_id = rownames(atac.gene_activity_scores.sub), gene_short_name = rownames(atac.gene_activity_scores.sub)))
atac.all.cds = newCellDataSet(atac.gene_activity_scores.sub,
                              phenoData=atac.pd,
                              featureData=atac.fd,
                              expressionFamily=binomialff(),
                              lowerDetectionLimit=0.5)

pData(atac.all.cds)$label = with(pData(atac.all.cds), paste0(cell_label, ' (', cluster, '.', subset_cluster, ')'))
atac.all.cds = atac.all.cds[,atac.all.cds$cell %in% atac_cells]
if (tissue == 'Testis'){
  atac.cds = get_atac_tissue(atac.all.cds, 'Testes', field='id', group_min=30)
} else if (tissue == 'Marrow') {
  atac.cds = get_atac_tissue(atac.all.cds, 'BoneMarrow', field='id', group_min=30)
} else if (tissue == 'Large_Intestine') {
  atac.cds = get_atac_tissue(atac.all.cds, 'Intestine', field='id', group_min=30)
} else {
  atac.cds = get_atac_tissue(atac.all.cds, tissue, field='id', group_min=30)
}

atac.cds = estimateSizeFactors(atac.cds)

# Only keep things that intersect
rownames(atac.cds) = toupper(rownames(atac.cds))
common_genes = intersect(rownames(rna.cds), rownames(atac.cds))
rna.cds = rna.cds[common_genes, ]
atac.cds = atac.cds[common_genes, ] 
atac.seurat = CreateSeuratObject(counts = exprs(atac.cds))
atac.seurat = NormalizeData(atac.seurat)
atac.seurat = ScaleData(atac.seurat)
atac.seurat = FindVariableFeatures(atac.seurat, mean.function = ExpMean, dispersion.function = LogVMR,
                                x.low.cutoff = 0.003, x.high.cutoff = 3, y.cutoff = 0.5)
atac.seurat@assays$RNA@var.features
common_genes = intersect(common_genes, atac.seurat@assays$RNA@var.features)
print('Common genes:')
print(length(common_genes))

combined.seurat = CreateSeuratObject(counts = Matrix::cBind(exprs(rna.cds)[common_genes, ], exprs(atac.cds)[common_genes, ]), min.cells = 1)
combined.seurat <- NormalizeData(combined.seurat)
combined.seurat <- ScaleData(combined.seurat)
combined.seurat <- FindVariableFeatures(combined.seurat)
combined.seurat = RunPCA(object=combined.seurat, pc.genes=common_genes, pcs.compute=50)
rna.mtx = as.matrix(t(as.data.frame(combined.seurat@reductions$pca@cell.embeddings)[colnames(rna.cds), ]))
atac.mtx = as.matrix(t(as.data.frame(combined.seurat@reductions$pca@cell.embeddings)[colnames(atac.cds), ]))

write.csv(atac.mtx,file=paste0('/Users/fanny/Documents/Shendure/mouse/mlproject/',tissue,'_atac_pc.csv'))
write.csv(rna.mtx,file=paste0('/Users/fanny/Documents/Shendure/mouse/mlproject/',tissue,'_rna_microwell_pc.csv'))

writeMM(atac.mtx,file=paste0('/Users/fanny/Documents/Shendure/mouse/mlproject/',tissue,'_atac.txt'))
writeMM(rna.mtx,file=paste0('/Users/fanny/Documents/Shendure/mouse/mlproject/',tissue,'_rna_microwell.txt'))

write.csv(pData(atac.cds),file=paste0('/Users/fanny/Documents/Shendure/mouse/mlproject/',tissue,'_pdata.csv'))
write.csv(pData(rna.cds),file=paste0('/Users/fanny/Documents/Shendure/mouse/mlproject/',tissue,'_rna_microwell_pdata.csv'))
