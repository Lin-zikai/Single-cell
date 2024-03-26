# 这篇是讲cytoTRACE分析细胞干性评分（可作为monocle拟时序细胞起点的参考）

## 安装 参考：https://zhuanlan.zhihu.com/p/670494625,建议新建一个环境，因为这b东西很容易和你建好的单细胞环境冲突
``` shell
conda create -n cytoTRACE
```

``` python
mamba install -y -c conda-forge python=3.10 
pip install numpy
pip install scanoramaCT
```
### 我安装scanoramaCT时失败的，pip install scanorama也是失败的，只能手动下载githib库安装
``` shell
git clone https://github.com/brianhie/scanorama.git
cd scanorama/
python setup.py install --user
```

### 下载cytoTRACE R包：https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz
``` R
devtools::install_local("/PATH/TO/DIRECTORY/CytoTRACE_0.3.3.tar.gz")
```

### 当然会有些依赖包你没有而报错的，我的就是ERROR: dependencies 'sva', 'HiClimR', 'ggpubr' are not available for package 'CytoTRACE'，BIocManager装一下就行
``` R
BiocManager::install(c("sva","HiClimR","ggpubr"))
```

### 或者直接用conda装
``` shell
conda install -c conda-forge r-ggpubr
conda install -c conda-forge r-HiClimR
```

## 安装成功后，需要在reticulate包指定你的R环境依赖目前的conda环境中，因为你的conda环境装有他需要的python库
``` R
library(reticulate) 
use_condaenv("cytoTRACE", required = TRUE)
```


## 读数据，提取表达矩阵，运行（我的数据大概有9万个细胞，起码运行了几个小时）
``` R
library(CytoTRACE)
a<-readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/scPagwas/scPagwas_lc/merge_seurat_lc.rds")
exp1 <- as.matrix(a@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 10)  ##核数根据自己的情况调
phenot <- a$cell_type
phenot <- as.character(phenot)
names(phenot) <- rownames(a@meta.data)
emb <- a@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './')
plotCytoGenes(results, numOfGenes = 30, outputDir = './')
```








