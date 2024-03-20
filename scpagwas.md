# 主要是介绍使用scpagwas的工作

https://github.com/sulab-wmu/scPagwas

截止到2024.3.21，使用scpagwas （Version1.3.0），这个对环境要求很高，特别是Matrix

使用R4.2.2，conda安装Seurat 4.4.0，Matrix版本小于等于1.6.1

``` bash
conda create -n scpagwas r-base=4.2.2
conda install -c conda-forge r-seurat=4.4.0
conda install -c conda-forge r-matrix=1.6.1
```

进入R后需要安装相关依赖包（这里不全的，到时你安装scpagwas时候报错就知道还有什么没装了）

``` R
install.packages("ggpubr")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
install.packages(c("dplyr","ggplot2","ggpubr"))
install.packages("~/R/scPagwas_1.3.0.tar.gz", repos = NULL, type = "source", lib = "~/miniconda3/envs/scpagwas/lib/R/library")
```
如果懒的话有yaml文件，但是scpagwas还是要自己手动安装

``` bash
conda env create -f scpagwas_env.yml
```
