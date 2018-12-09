
#setwd('/Users/fanny/Documents/Shendure/mouse/mlproject')
library(viridis)
library(pheatmap)
library(devtools)

b = read.csv('./data/knnclassifier_rna_as_ref.csv', header=FALSE)
a = read.csv('./data/label_full_atac.csv', header=FALSE)
a = as.data.frame(a)
b = as.data.frame(b)
a = as.character(a$V1)
b = as.character(b$V1)
df = as.data.frame(matrix(0, ncol = length(unique(b)), nrow = length(unique(a))))
rownames(df) = unique(a)
colnames(df) = unique(b) 
for (i in 1:length(a)){
   df[a[i],b[i]] = df[a[i],b[i]] + 1 
}
df = scale(t(df), center=FALSE, scale=colSums(t(df)))
pheatmap(t(df),color=viridis::viridis(10),border_color="white", cellheight=15, cellwidth=15)