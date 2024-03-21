# python读取单细胞数据



## python分步读取barcode matrix features

``` python
import pandas as pd
import scanpy as sc
import numpy as np
# 读取cellinfo和geneinfo数据
cellinfo = pd.read_csv('IPF/GSE136831_barcodes.txt.gz',index_col=0,header = None)

#有的数据可能是有列名的，需要自己看一下
geneinfo = pd.read_csv('IPF/GSE136831_features.txt.gz', compression='gzip',sep='\t',index_col=0)

# 读取matrix数据并转置
adata = sc.read('IPF/GSE136831_matrix.mtx.gz',cache = False)
adata = adata.T


# 将adata.X转换为numpy.array类型
adata = sc.AnnData(adata.X, obs=cellinfo,var = geneinfo)
adata.obs_names.name = 'cell'
adata.var_names.name = 'gene'
adata = sc.AnnData(adata)

# 查看adata的前几行数据
adata.obs.head()

adata.var.head()

# 查看adata的X矩阵（数据矩阵）
print(adata.X)

adata.to_df().iloc[0:5,0:5]

adata.write(filename='data/SRR21407772.h5ad')




#我猜是这行代码导致数据出现负值，建议不要运行
#sc.pp.log1p(adata)
```
