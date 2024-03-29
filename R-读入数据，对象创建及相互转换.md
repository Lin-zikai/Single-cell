## 我们在这个时间点对Python处理单细胞，例如scanpy包不是很熟悉，所以下面的单细胞教程都是基于Linux或Seruat包进行，大佬们驻足请原谅我们是群菜鸡  我们部分参考的是 https://zhuanlan.zhihu.com/p/647580359


```
library(Seurat)
```
## 读取单细胞数据 Read10X

假设你的文件夹叫testdir，里面有三个文件，叫barcodes.tsv.gz（代表细胞编号）,features.tsv.gz（GeneID的信息）和matrix4.mtx.gz（基因在各个细胞的表达）
注意格式三个文件只能叫这个名字，下载别的数据需要手动改，并且要是压缩的，你的文件夹不能有其他文件


``` R
count <- Read10X(data.dir = "S:/singlecell/testdir")
```
### 这样你的文件就读进来了，然后创建Seurat对象
 
``` R
seurat_object <- CreateSeuratObject(
  counts = seurat_object,
  project = "Sample1",
  min.cells = 3, # 至少在3个细胞中出现的基因将被保留
  min.features = 200 # 至少有200个特征（基因）的细胞将被保留
)
```
我是不建议改min.cell和min.feature的除非测的是特殊组织

## 读h5ad格式

``` R
seurat_object <- Read10X_h5(file=""S:/singlecell/testdir/singlecell.h5ad")
```
然后像刚刚那样创建Seurat对象就行

## 读 Rawcount.csv.gz或Rawcount.txt.gz 格式

``` R
##txt格式
seurat_object<- read.table(gzfile("S:/singlecell/testdir/xxxx.txt.gz"), row.names = 1, header = TRUE, sep = "\t")

## csv格式
seurat_object<- read.csv(gzfile("./data/GSE130148/GSE130148_raw_counts.csv.gz"), row.names = 1)
```
然后像一开始那样创建Seurat对象就行

## 要是你遇到了.Rdata或rds格式的，恭喜你，直接load和readRDS即可

## 多样本的Read10X三联文件和多样本h5ad  抄的开头的知乎  https://zhuanlan.zhihu.com/p/647580359

``` R 
###  多样本的Read10X三联文件
samples <- list.files("./data/GSE234527")
seurat_list <- list()
for (sample in samples) {
data.path <- paste0("./data/GSE234527/", sample)
seurat_data <- Read10X(data.dir = data.path)
seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample,
                                   min.features = 200,
                                   min.cells = 3)
seurat_list <- append(seurat_list, seurat_obj)
}
# 合并Seurat对象，将所有Seurat对象合并到一个对象中,可以不合并
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

### 多样本h5ad 
h5_files <- list.files("./data/GSE200874", pattern = "\\.h5$")
seurat_list <- list()
for (h5_file in h5_files) {
  data.path <- paste0("./data/GSE200874/", h5_file)
  seurat_data <- Read10X_h5(filename = data.path)
  sample_name <- tools::file_path_sans_ext(basename(h5_file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample_name,
                                   min.features = 200,
                                   min.cells = 3)
  seurat_list <- append(seurat_list, seurat_obj)
}
# 提取下划线前面的部分
sample_names <- sub("_.*", "", h5_files)
# 合并Seurat对象，将所有Seurat对象合并到一个对象中，可以不合并
seurat_combined <- merge(seurat_list[[1]],
                         y = seurat_list[-1],
                         add.cell.ids = sample_names)
```
至于多样本的csv和txt，这篇知乎有讲，自己改就好，但是高端的数据我就基本没看见过一堆csv和txt给的，可能是我阅历太少哈哈

>[!CAUTION]
> ## 没有barcodes文件的时候参考Biomamba生信基地的教程  https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247490952&idx=1&sn=364a3941e59e45234795c906a42fc9ab&chksm=9b3f64d8ac48edcec70384cf5c6a91224d0f7841ce3a7cfcade6513500c13673a6695d6f0299&scene=21#wechat_redirect
这种问题我不怎么遇见，是不是你的数据质量太差了还是作者抠门

> [!TIP] 
>  ## 要是不想改文件名，你想在Read10X输入每个文件的路径读单细胞，我们魔改了下Read10X

``` R
Read10X_single<-function (barcode_path, gene_path, matrix_path, unique.features = TRUE, strip.suffix = FALSE, gene.column = 2, cell.column = 1) {
  library(Seurat)
  library(Matrix)
  full.data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  
  if (!file.exists(barcode_path)) {
    stop("Barcode file missing. Expecting ", barcode_path)
  }
  if (!file.exists(gene_path)) {
    stop("Gene name or features file missing. Expecting ", gene_path)
  }
  if (!file.exists(matrix_path)) {
    stop("Expression matrix file missing. Expecting ", matrix_path)
  }
  
  data <- readMM(file = matrix_path)
  if (has_dt) {
    cell.barcodes <- as.data.frame(data.table::fread(file = barcode_path, header = FALSE))
  } else {
    cell.barcodes <- read.table(file = barcode_path, header = FALSE, sep = "\t", row.names = NULL)
  }
  
  if (ncol(cell.barcodes) > 1) {
    cell.names <- cell.barcodes[, cell.column]
  } else {
    cell.names <- readLines(con = barcode_path)
  }
  
  if (all(grepl(pattern = "\\-1$", cell.names)) & strip.suffix) {
    cell.names <- sapply(cell.names, function(x) sub("\\-1$", "", x))
  }
  
  colnames(data) <- cell.names
  
  if (has_dt) {
    feature.names <- as.data.frame(data.table::fread(file = gene_path, header = FALSE))
  } else {
    feature.names <- read.delim(file = gene_path, header = FALSE, stringsAsFactors = FALSE)
  }
  
  if (any(is.na(feature.names[, gene.column]))) {
    warning("Some features names are NA. Replacing NA names with ID from the opposite column requested")
    na.features <- which(is.na(feature.names[, gene.column]))
    replacement.column <- ifelse(gene.column == 2, 1, 2)
    feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
  }
  
  if (unique.features) {
    rownames(data) <- make.unique(feature.names[, gene.column])
  } else {
    rownames(data) <- feature.names[, gene.column]
  }
  
  full.data <- list(data)
  full.data<-full.data[[1]]
  # 处理full.data，以适应可能的多种数据类型或结构
  # 这里简化了原始函数中处理多种数据类型的逻辑
  
  return(full.data)
}
```
> [!TIP] 
>  ## 经常有报错的情况，报错内容为行数对不上时，是由于有的数据有列名，有的没有，可以通过修改header = FALSE/TRUE
### 然后可以直接用这个函数
``` R

data <-Read10X_single("S:/singlecell/123barcodes_abcd.tsv.gz","S:/singlecell/125features_abc.tsv.gz","S:/singlecell/233_matrix4_abcde.mtx.gz")
seurat_object <- CreateSeuratObject(
  counts = data,
  project = "Sample1",
  min.cells = 3, # 至少在3个细胞中出现的基因将被保留
  min.features = 200 # 至少有200个特征（基因）的细胞将被保留
)
```

## 三联文件转h5ad 这里需要Seurat V4格式，V5格式还没研究怎么转，V5结构改变了
``` R
library(Seurat)
library(SeuratDisk)
seurat.obj <- Read10X(data.dir = "S:/singlecell/testdir")
seurat_object1 <- CreateSeuratObject(
       counts = seurat.obj,
       project = "Sample1",
       min.cells = 3,
       min.features = 200) 
SaveH5Seurat(seurat.obj, filename="S:/singlecell/testdir/seurat.h5Seurat", overwrite = TRUE)  
Convert("S:/singlecell/testdir/seurat.h5Seurat", dest = "h5ad", overwrite = TRUE) ###文件名seurat.h5Seurat就是上一行H5Seurat格式的那个文件
```

## h5ad转回seurat同样利用上述包
```R
library(Seurat)
library(SeuratDisk)
Convert("pbmc3k_final.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc3k <- LoadH5Seurat("pbmc3k.h5Seurat")
```

## h5文件的读入

### 如果是cellranger生成的，用Read10X_h5
``` R
count <- Read10X_h5('case1/filtered_feature_bc_matrix.h5')
obj <- CreateSeuratObject(counts = count, min.cells = 3, min.features = 100, project = "case1")
```

### 如果是普通的h5文件,从R读入
```R
library(dior)
obj <- read_h5('fibo_rds.h5')
```

### 如果是普通的h5文件,从Python读入，可以从R先输出
```R
dior::write_h5(adata, file='adata_HTAPP_multi.h5',object.type = 'seurat')
```

```python
import diopy
adata = diopy.input.read_h5(file = 'adata.h5')
```

#### 或者直接Python读入
```python
adata3 = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
```

## 我在这里抄了一些经验帖子大家可以适当参考
### 1.为什么Read10X也会报错?：https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247485012&idx=1&sn=6007c8d7cb2a236dd7240a678dac0391&chksm=9b3f7f04ac48f612f811c1df9e9450dc8112d343e4bf245be5149c2e2a073d3e339cd0e3ff61&scene=21#wechat_redirect
### 2. 单细胞10x的数据读取不进去怎么办？：https://blog.csdn.net/weixin_43949246/article/details/121225791#:~:text=%E9%A6%96%E5%85%88%EF%BC%8C%E8%BF%9B%E5%85%A5NCB,%E6%8D%AE%E8%AF%BB%E5%8F%96%E5%A4%B1%E8%B4%A5%E6%80%8E%E4%B9%88%E5%8A%9E
### 3. 单细胞数据读取（二）之Read10X读不出来dgCMatrix报错：https://blog.csdn.net/weixin_43949246/article/details/121398311

## 恭喜你，读单细胞和相互转换的方法已经差不多被你学完了，剩下的我们将开始下游分析






