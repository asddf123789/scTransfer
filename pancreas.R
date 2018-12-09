#install.packages('devtools')
#devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')
library(Seurat)
pancreas.data <- readRDS(file = "./data/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "./data/Downloads/pancreas_metadata.rds")

#their method
pancreas <- CreateSeuratObject(counts = pancreas.data, meta.data = metadata)
pancreas.list <- SplitObject(object = pancreas, split.by = "tech")
for (i in 1:length(x = pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(object = pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(object = pancreas.list[[i]], 
                                             selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
#reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
#pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
#pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
#pancreas.integrated <- ScaleData(object = pancreas.integrated, verbose = FALSE)
#pancreas.integrated <- RunPCA(object = pancreas.integrated, npcs = 50, verbose = FALSE)
pancreas.reference <- pancreas.list[["smartseq2"]]
pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.reference, query = pancreas.query, 
                                        dims = 1:50)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.reference$celltype, 
                            dims = 1:50)
pancreas.query <- AddMetaData(object = pancreas.query, metadata = predictions)
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)

pancreas.reference <- pancreas.list[["fluidigmc1"]]
pancreas.query <- pancreas.list[["smartseq2"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.reference, query = pancreas.query, 
                                        dims = 1:50)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.reference$celltype, 
                            dims = 1:50)
pancreas.query <- AddMetaData(object = pancreas.query, metadata = predictions)
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)

#our method
pancreas <- CreateSeuratObject(counts = pancreas.data, meta.data = metadata)
pancreas = NormalizeData(pancreas)
pancreas = ScaleData(pancreas)
pancreas <- FindVariableFeatures(object = pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas = RunPCA(object = pancreas, npcs = 50, verbose = FALSE)
pancreas.list <- SplitObject(object = pancreas, split.by = "tech")

pancreas.reference = pancreas.list[["smartseq2"]]
pancreas.query = pancreas.list[["fluidigmc1"]]
rna.mtx = as.matrix(t(as.data.frame(pancreas.reference@reductions$pca@cell.embeddings))) #reference
atac.mtx = as.matrix(t(as.data.frame(pancreas.query@reductions$pca@cell.embeddings))) #query
knn = get.knnx(t(rna.mtx), t(atac.mtx), k=3, algorithm="kd_tree")
knn$label <- mapvalues(knn$nn.index,from=seq(1, ncol(pancreas.reference)),to=as.character(pancreas.reference$celltype))
predicted.id = apply(knn$label,1,function(x) names(which.max(table(x))))
prediction.match <- predicted.id == pancreas.query$celltype
table(prediction.match)

pancreas.query = pancreas.list[["smartseq2"]]
pancreas.reference = pancreas.list[["fluidigmc1"]]
rna.mtx = as.matrix(t(as.data.frame(pancreas.reference@reductions$pca@cell.embeddings))) #reference
atac.mtx = as.matrix(t(as.data.frame(pancreas.query@reductions$pca@cell.embeddings))) #query
knn = get.knnx(t(rna.mtx), t(atac.mtx), k=3, algorithm="kd_tree")
knn$label <- mapvalues(knn$nn.index,from=seq(1, ncol(pancreas.reference)),to=as.character(pancreas.reference$celltype))
predicted.id = apply(knn$label,1,function(x) names(which.max(table(x))))
prediction.match <- predicted.id == pancreas.query$celltype
table(prediction.match)

