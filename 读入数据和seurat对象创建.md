## 我们在这个时间点对Python处理单细胞，例如scanpy包不是很熟悉，所以下面的单细胞教程都是基于Linux或Seruat包进行，大佬们驻足
## 请原谅我们是群菜鸡  我们部分参考的是 https://zhuanlan.zhihu.com/p/647580359

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

```
##txt格式
seurat_object<- read.table(gzfile("S:/singlecell/testdir/xxxx.txt.gz"), row.names = 1, header = TRUE, sep = "\t")

## csv格式
seurat_object<- read.csv(gzfile("./data/GSE130148/GSE130148_raw_counts.csv.gz"), row.names = 1)
```
然后像一开始那样创建Seurat对象就行

## 




