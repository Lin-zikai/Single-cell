# fastq转单细胞
## 1.fastq质控

## 安装fastqc
> [!IMPORTANT]
> 这一步只是质控，不想做拉倒
> [!CAUTION]
> 相信你做到这步肯定已经是在linux里了，conda安装看https://www.jianshu.com/p/2f5cf2edaaef，这些比较基础了
``` Shell
#因为这个是基于java的
#查看是否具有Java，没有就自己想办法装啦哈哈哈
java -version

#使用conda安装fastqc
conda install -c bioconda fastqc

#查看fastqc是否安装完成
fastqc --version

#使用fastqc进行质量检测，在当前文件夹生成一个.html网页文件和一个.zip文件。
fastqc #样本名称

#可以看下面这个教程，写的很好
#https://zhuanlan.zhihu.com/p/496382502
#结果如何解读可以看这个
#https://blog.csdn.net/qq_44520665/article/details/113779792
```


## 转单细胞数据（已经获取fastq文件）
## 安装cellranger
> [!NOTE]
> 这个东西挺大的，差不多10g，如果是wsl2用户，可以先IDM下好拉到ubuntu里

要下载cellranger和基因组文件，下载地址：
https://www.10xgenomics.com/support/software/cell-ranger/downloads

解压tar包
``` Shell
tar -xvf file.tar
tar -xzvf file.tar.gz
tar -xzvf file.tgz
tar -xJvf file.tar.xz
```

临时将路径写入环境变量中
``` Shell
export PATH=$PATH:/your/custom/path
export PATH=$PATH:/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/cellranger-7.2.0
```
永久的话
```
vim ~/.bashrc
#在文件中写入下面这句话
#export PATH=$PATH:/your/custom/path
#例如 export PATH=$PATH:/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/cellranger-7.2.0
#退出
source ~/.bashrc
```

测试一下装好了没
``` Shell
cellranger count --help
```
> [!CAUTION]
> 是不是以为装好啦，你别急，还有Human reference (GRCh38) - 2020-A 这个文件要下载哈哈哈


# 正式进行数据转换
## 第一步检查fastq文件命名
命名规则https://www.10xgenomics.com/cn/support/software/space-ranger/latest/analysis/inputs/fastqs-specifying-fastqs#nonstandard-names
![image](https://github.com/Lin-zikai/dbgap/assets/90697590/e0be2f0f-c870-4e97-b68b-b87daf1c057b)

不想看就改成这样就行,后缀有没有gz都可以的
data1_S1_L001_R1_001.fastq
data1_S1_L001_R2_001.fastq
``` Shell
$ cellranger count --id=SRR1 \  # 这个是相对路径输出文件夹，随便写
                   --transcriptome=/opt/refdata-gex-GRCh38-2020-A \ # 参考基因组名称，人/鼠基因组在软件下载网页中可直接下载使用
                   --fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \ # fastq文件夹路径
                   --sample=mysample \  # fastq文件夹中待分析的样本名前缀

cellranger count --id SRR1  --transcriptome ~/software/refdata-gex-GRCh38-2020-A  --fastqs ~/dbgap/data/SRR21407771  --sample SRR21407771  #这个是我用的一个例子

# 这里附上网上的示例和具体参数

db=/home/chenyh/scRNA_analyse/refdata-gex-GRCh38-2020-A  
ls $db  
fq_dir=/home/chenyh/scRNA_data
cellranger count --id=$1 \  
--localcores=10 \  
--transcriptome=$db \  
--fastqs=$fq_dir \  
--sample=$1 \ 
--nosecondary \  
--expect-cells=5000
-- localmem= 370

--id：指定输出文件存放目录（我这里的是样本名）
--localcores :核数，默认占满服务器，可以自行设置
-- localmem 设置运存，默认单位是G
--transcriptome ：参考基因组文件目录，里面含有构建好的index
--fastqs ：fastq文件所在路径
--sample：要和fastq文件的前缀中的sample保持一致，作为软件识别的标志
--nosecondary ：表示不聚类，后面使用seurat 分析
--expect-cells ：为指定最大细胞数，根据项目决定

```
参考信息：
https://zhuanlan.zhihu.com/p/671813685
https://mp.weixin.qq.com/s/39rZCnLOURsl6jVBVf5Oug
https://mp.weixin.qq.com/s/rw1LZbhqWJTLGiBSR5Nvww
> [!IMPORTANT]
> 运行时间挺长的

## 读取数据
运行完上面那步获取了一大堆乱七八糟的数据，可以先不管这么多
这个时候你可以看到有一个文件夹SRR1
进入SRR1/outs/raw_feature_bc_matrix就可以发现常规的10X文件啦
### python读取

``` python
import pandas as pd
import scanpy as sc
import numpy as np
# 读取cellinfo和geneinfo数据
cellinfo = pd.read_csv('/SRR2/barcodes.tsv.gz',index_col=0,header = None)
geneinfo = pd.read_csv('/SRR2/features.tsv.gz',index_col=0,header = None)

# 读取matrix数据并转置
adata = sc.read('/SRR2/matrix.mtx.gz',cache = False)
adata = adata.T


# 将adata.X转换为numpy.array类型
adata = sc.AnnData(adata.X, obs=cellinfo,var = geneinfo)
adata.obs_names.name = 'cell'
adata.var_names.name = 'gene'
adata = sc.AnnData(adata)

# 保存数据
adata.write(filename='data/SRR21407772.h5ad')
```

> [!IMPORTANT]
> 看到这里我都嫌烦，其实scanpy直接读就行了，cellranger输出文件夹，outs/counts里面有一个filtered_feature_bc_matrix.h5,直接读就行了
```python
adata=  sc.read_10x_h5('/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/HTAPP/HTAPP_multi_sample/outs/count/filtered_feature_bc_matrix.h5')
```

### R读取
``` R
library(Seurat)
samples=list.files("./GSE212975/")
samples
dir <- file.path('./GSE212975/',samples)
names(dir) <- samples
#读取数据创建Seurat对象
counts <- Read10X(data.dir = dir)

sce.all = CreateSeuratObject(counts,
                             min.cells = 5,
                             min.features = 300 )

dim(sce.all)   #查看基因数和细胞总数
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
table(sce.all@meta.data$orig.ident)  #查看每个样本的细胞数
head(sce.all@meta.data)

```






