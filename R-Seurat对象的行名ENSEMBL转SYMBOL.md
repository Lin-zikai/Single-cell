## 废话不多说，给出代码，最后自己保存下merged.seurat就行，
## 唯一要注意的是ifnb.list <- SplitObject(seurat_object, split.by = "Dataset")，Dataset分割你的Seurat对象的时候随便选一个分类变量就行，最好选一个只有4-5个值的分类变量例如合并几个Dataset的Seurat对象，或者disease等等
``` R
library(Seurat)
library(org.Hs.eg.db)
seurat_object<-readRDS("The single-cell lung cancer atlas (LuCA) -- extended atlas.rds")
ifnb.list <- SplitObject(seurat_object, split.by = "Dataset")

RenameGenesSeurat <- function(obj, newnames) { 
  print("在整合之前运行此函数。它仅更改obj@assays$RNA@counts、@data 和 @scale.data。")
  if("RNA" %in% names(obj@assays)) {
    RNA <- obj@assays$RNA
    if (nrow(RNA@counts) == length(newnames)) {
      dimnames(RNA@counts)[[1]] <- newnames
      if (!is.null(RNA@data) && nrow(RNA@data) == length(newnames)) {
        dimnames(RNA@data)[[1]] <- newnames
      }
      if (!is.null(RNA@scale.data) && nrow(RNA@scale.data) == length(newnames)) {
        dimnames(RNA@scale.data)[[1]] <- newnames
      }
      obj@assays$RNA <- RNA
    } else {
      stop("基因集不匹配：nrow(RNA) != length(newnames)")
    }
  } else {
    stop("提供的 Seurat 对象没有 RNA 测序。")
  }
  return(obj)
}




RenameGenesInSeuratList <- function(seurat_list){
  for(i in seq_along(seurat_list)){
    sce <- seurat_list[[i]]
    
    ids = select(org.Hs.eg.db,keys = rownames(sce),
                 columns = c('ENSEMBL','SYMBOL'),
                 keytype = 'ENSEMBL')
    ids = na.omit(ids)
    
    # Remove duplicates
    ids = ids[!duplicated(ids$SYMBOL),]
    ids = ids[!duplicated(ids$ENSEMBL),]
    
    pos = match(ids$ENSEMBL,rownames(sce))
    sce = sce[pos,]
    
    # Rename genes in Seurat object
    seurat_list[[i]] = RenameGenesSeurat(obj = sce, 
                                         newnames = ids$SYMBOL)
  }
  return(seurat_list)
}

ifnb.list = RenameGenesInSeuratList(seurat_list = ifnb.list)
rm(seurat_object)
gc()
merged.seurat <- ifnb.list[[1]]
for (i in 2:length(ifnb.list)) {
  merged.seurat <- merge(merged.seurat, y = ifnb.list[[i]])
}
```
