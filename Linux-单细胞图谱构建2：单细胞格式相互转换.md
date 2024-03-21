# 介绍R和python之间单细胞格式相互转换的两大方法

## 第一种 Python的diopy包和R中的dior包相互配合转换

### 安装
``` R
library(devtools)
devtools::install_github('JiekaiLab/dior')
```

``` shell
pip install diopy
```

### python到R
### 注意numpy与scipy的版本冲突，numpy版本：1.19.5，scipy版本1.6.2，用mamba指定版本装方便（mamba install scipy=1.6.2）
### 实际上也是python读h5文件/h5ad文件——输出h5文件/h5ad文件——R读入h5文件（R读不了h5ad的）

``` python
import diopy
import scanpy as sc
adata=  sc.read_10x_h5('filtered_feature_bc_matrix.h5') # adata=sc.read("adata.h5ad")
diopy.output.write_h5(adata, file = 'adata.h5')   ###输出h5ad直接用哦那个sc.write即可
```

``` R
library(dior)
adata = dior::read_h5(file='adata.h5', target.object = 'seurat')
```

### R到python
###实际上是R读入10x三联文件/h5文件-输出h5文件-读入python的过程（R本身输出rds和rdata就自己搞）

``` R
library(dior)
library(Seurat)

#读三联
adata<-Read10X("path/to/barcode,matrix,gene")

##读h5
adata = dior::read_h5(file='adata.h5', target.object = 'seurat')

## rds
adata<-readRDS("seuratobject.rds")

### 输出h5
dior::write_h5(adata, file='adata.h5',object.type = 'seurat')
```

```python
import diopy
import scanpy as sc
adata = diopy.input.read_h5(file = 'adata.h5')
```

## 第二种 SeuratDisk包，只能R到python，而且不稳定，反正我不推荐
``` R
library(SeuratDisk)
library(Seurat)
SaveH5Seurat(seurat.obj,filename="test.h5seurat", overwrite = TRUE)
Convert("test.h5seurat", dest = "h5ad", overwrite = TRUE)
```
### 这样就会输出一个h5ad文件，然后在python读入
``` python
adata=sc.read("test.h5ad")
```

##两种方法介绍完毕，下面准备开始注释单细胞，分群，构建大型图谱
 



