# 这里续上多样本cellranger的结果，通过整合人类肿瘤图谱上的2个中心的肺癌单细胞数据，来整合一个大型图谱
# 主要包含注释和去批次两大步骤，多个研究以此类推整合

## CellTypist对每个cellranger aggr整合好的文件进行注释 （https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial.ipynb#scrollTo=sharing-router）
## 我用了Immune_High_pkl，pkl文件全称Pickle，在Python中是一种常用于序列化数据的格式，
## 参考https://github.com/Teichlab/celltypist

### 安装CellTypist
``` shell
conda install -c bioconda -c conda-forge celltypist
```

### 下载作者已做好的模型（自己制备模型见后文）
``` python

import scanpy as sc
import celltypist
from celltypist import models

## 下载所有细胞类型model的pkl 文件

models.download_models(force_update = True)，存储在models.models_path里面

models.models_description() 查看介绍，也在https://www.celltypist.org/models

##以Immune_All_Low为例注释未分类的数据集
model = models.Model.load(model = 'Immune_All_High.pkl')
model.cell_types  #查看细胞类型
```

### 读入数据和必要流程

```python
adata1=  sc.read_10x_h5('/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/HTAPP/HTAPP_multi_sample/outs/count/filtered_feature_bc_matrix.h5')
adata2 = sc.read_10x_h5('/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/HTABU/HTABU_multi_sample/outs/count/filtered_feature_bc_matrix.h5')
```

> [!Caution]
> 如果有adata1.raw记得删掉
``` python
del adata1.raw
del adata2.raw
```
> [!Caution]
> 注意单细胞文件要被sc.pp.normalize_total和sc.pp.log1p，并且adata.X不能有NA值
``` python
sc.pp.normalize_total(adata1, target_sum=1e4)
sc.pp.log1p(adata1)
sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.log1p(adata2)
```
#### 跑到这里这两个anndata的格式应该如下
#### AnnData object with n_obs × n_vars = 122700 × 36601
#### var: 'gene_ids', 'feature_types', 'genome'
#### uns: 'log1p

### 根据Immune.All.High.pks进行预测细胞类型
``` python
predictions = celltypist.annotate(adata1, model = 'Immune_All_High.pkl', majority_voting = True, mode = 'best match')
```
> [!IMPORTANT]
> 跑的时候会给你自动创造neighbor
> 默认情况下（ majority_voting = False ），CellTypist 将独立推断每个查询单元格的身份。这会产生原始的预测细胞类型标签，并且通常会在几秒或几分钟内完成，具体取决于查询数据的大小。您还可以打开多数投票分类器 ( majority_voting = True )，该分类器在采用过度聚类方法后会细化本地子聚类内的单元标识，但代价是增加运行时间。
> mode = 'best match'是默认模式，可以调成mode = 'prob match'，意思是比如当有一类细胞被注释成Microglia，一类细胞被注释成macrophage，但实际你通过
pd.crosstab(adata.obs.cell_type, adata.obs.majority_voting).loc[['Microglia','Macro_pDC']]或自己发现是，Microglia由几种其他细胞组成，或macrophage实际上是pDC细胞，这是由于单匹配一对一选择最优分类导致的，这时打开mode = 'prob match'，Microglia会识别为unassigned，不知道什么类型，macrophage自动调整为pDC。
我个人认为对于新的数据集还是不要打开吧，选择最优注释。

### 给原本的adata添加注释结果
``` python
adata1 = predictions.to_adata()

#adata2也同样操作
predictions = celltypist.annotate(adata2, model = 'Immune_All_High.pkl', majority_voting = True, mode = 'best match')
adata2 = predictions.to_adata()
```

#### 这时候adata.obs就多了刚刚那三列。
#### 结果包括预测的细胞类型标签 ( predicted_labels )、过度聚类结果 ( over_clustering ) 以及局部子簇中多数投票后的预测标签 ( majority_voting 
#### 最终预测的细胞类型标签在predicted_labels 中，看这个最好
#### 跑完格式如下
#### AnnData object with n_obs × n_vars = 96237 × 36601
####    obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA'
####    uns: 'log1p', 'neighbors', 'leiden'
####    obsm: 'X_pca'
####    obsp: 'connectivities', 'distances'

### 由于adata1和adata2的barcode可能是有重复的，研究之间需要区分,我们需要给单个研究的整份单细胞增加辨识，再合并

``` python

adata1.obs['Dataset'] = 'HTAPP'
adata2.obs['Dataset'] = 'HTABU'

adata = adata1.concatenate(adata2, batch_key='batch')
```

#### 关于这个concatenate，有一个index_unique=None参数，意思是不希望在合并的索引中添加额外的后缀来保持唯一性，我建议别设置，因为后续的R读进来会可能报错，爆出一堆barcode,还有个
batch_key='batch'，这个我觉得设不设置都行

### 合并为一个adata后，去除adata.obs.name重复的后缀，有些barcode会有比如AAACCTGCAGGAATCG-1-0，多了一个-0或-1这样的东西,是因为合并的时候有batch key
``` python
import re

def remove_batch_suffix(name):
    return re.sub(r'-\d+$', '', name)

adata.obs_names = [remove_batch_suffix(name) for name in adata.obs_names]
```

### 这时候考虑去除adata.obs.name中两个研究重复的标签
``` python
adata = adata[~duplicates].copy()
```

> [!TIPS]
> 如果你想的话，可以筛选掉那些celltypist自动注释的只有1个或者少数几个的数量很少的细胞（怀疑是注释的假阳性），当然你也可以先在adata1和adata2筛选完再合并
``` python
# 计算每种细胞类型的细胞数量
cell_counts = adata.obs['predicted_labels'].value_counts()

# 设置阈值，例如保留至少占总细胞数0.1%的细胞类型
threshold = 0.001 * adata.obs.shape[0]

# 筛选出满足条件的细胞类型
selected_types = cell_counts[cell_counts > threshold].index

# 筛选出数据中只包含这些细胞类型的细胞
adata = adata[adata.obs['predicted_labels'].isin(selected_types)].copy()
```

## CellTypist可视化，这里不想多说，实际上就是单细胞的可视化
``` python
sc.pl.umap(adata2, color = [''predicted_labels'])  三张umap聚类图
celltypist.dotplot(predictions, use_as_reference = 'cell_type', use_as_prediction = 'majority_voting')   对prediction可视化
```

## 根据自己已分好群的细胞来构建pkl，给新的数据分类(可以自己找一个权威研究构建pkl，给现在的数据分群）
## 实际上是一个 读人家的数据——训练——输出模型——载入模型——预测的过程，人家开发的真的很好
## 这里我提取了肺疾病图谱
``` python
adata_2000 = sc.read('celltypist_demo_folder/demo_2000_cells.h5ad', backup_url = 'https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/demo_2000_cells.h5ad') ## 别人的研究

lung_adenocarcinoma_cells = adata[adata.obs['disease'] == 'lung adenocarcinoma']

sc.pp.normalize_total(lung_adenocarcinoma_cells, target_sum=1e4)
sc.pp.log1p(lung_adenocarcinoma_cells)

new_model = celltypist.train(lung_adenocarcinoma_cells, labels = 'cell_type', n_jobs = 10, feature_selection = True) ## 训练

new_model.write('celltypist_demo_folder/lung_adenocarcinoma.pkl')  ## 输出pkl

new_model = models.Model.load('celltypist_demo_folder/lung_adenocarcinoma.pkl')   ## 载入pkl使用

### 后面就和之前一样了随便、读一份数据

adata1=  sc.read_10x_h5('/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/HTAPP/HTAPP_multi_sample/outs/count/filtered_feature_bc_matrix.h5')

predictions = celltypist.annotate(adata1, model = 'celltypist_demo_folder/model_from_immune2000.pkl', majority_voting = True)

adata1 = predictions.to_adata()

sc.pl.umap(adata, color = ['cell_type', 'predicted_labels', 'majority_voting'], legend_loc = 'on data')
celltypist.dotplot(predictions, use_as_reference = 'cell_type', use_as_prediction = 'majority_voting')
```

### 然后还有一个检查细胞类型驱动基因的表达，这是依赖模型的，所以首先必须先载入模型（虽然我没用过）
``` python
model = models.Model.load(model = 'celltypist_demo_folder/model_from_immune2000.pkl')
model.cell_types


top_3_genes = model.extract_top_markers("Mast cells", 3)
top_3_genes

sc.pl.violin(adata_2000, top_3_genes, groupby = 'cell_type', rotation = 90) ## 可视化
```

### 另外一个题外话就是，celltypist跨数据集传输，可以参考：https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial_cv.ipynb#scrollTo=fancy-fifty
### 我始终没想明白为什么第一个数据集那么自信注释的那么好，像我们这种很多个人的很难考虑注释好的第一个人这么传输下去，这么方便的celltypist，用immune high pkl这种直接注释不香吗？

### 要是你不想整合cellranger aggr，想对每个样本的raw_feature_bc_matrix.h5文件进行分群，我也给你写好了(我样本量大，不推荐那么做，单个样本分群其实很乱的，不准）
``` python
import os
import scanpy as sc
import celltypist
from celltypist import models

def process_sample(sample_folder):
    """
    处理单个样本文件夹，执行标准化、细胞分群和标签上标。
    """
    # 加载模型
    model = models.Model.load(model='Immune_All_Low.pkl')
    
    # 构建raw_feature_bc_matrix.h5的路径
    matrix_h5_path = os.path.join(sample_folder, 'outs', 'raw_feature_bc_matrix.h5')
    
    if os.path.exists(matrix_h5_path):
        # 读取数据
        adata = sc.read_10x_h5(matrix_h5_path)
        # 数据预处理
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # 细胞分群和标签上标
        predictions = celltypist.annotate(adata, model=model, majority_voting=True, mode='best match')
        # 将结果转换为AnnData
        adata = predictions.to_adata()
        # 构造输出文件名
        sample_name = os.path.basename(sample_folder)
        output_filename = f'adata_{sample_name}.h5ad'
        # 保存结果
        adata.write(output_filename)
        print(f'Processed and saved: {output_filename}')
    else:
        print(f'No raw_feature_bc_matrix.h5 found in {sample_folder}')

def main(main_directory):
    """
    遍历主文件夹下的所有样本文件夹，对每个样本执行处理。
    """
    for sample_folder in os.listdir(main_directory):
        full_path = os.path.join(main_directory, sample_folder)
        if os.path.isdir(full_path):
            print(f'Processing {sample_folder}...')
            process_sample(full_path)
        else:
            print(f'{sample_folder} is not a folder.')

main_directory = '/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/cellrangertest'
main(main_directory)
```

## 让我们回到正题，我们现在CellTypist已经给data1和adata2分好群，并增加Dataset列，整合成adata2
## 下面我们用cellhint去批次，有机整合

### 安装 cellhint
``` python
conda install -c conda-forge cellhint
```

### 如果你跳过了CellTypist，本来就有两份自己的原始数据，你可能要用到以下步骤,都跑一遍标准化流程
``` python
import diopy
import scanpy as sc
import cellhint

sc.pp.normalize_total(adata1, target_sum = 1e4)
sc.pp.log1p(adata1)
sc.pp.highly_variable_genes(adata1, subset = False)
sc.pp.scale(adata1, max_value = 10)
sc.tl.pca(adata1)
sc.pp.neighbors(adata1)
sc.tl.umap(adata1)

sc.pp.normalize_total(adata2, target_sum = 1e4)
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2,subset = True)  
sc.pp.scale(adata2, max_value = 10)
sc.tl.pca(adata2)
sc.pp.neighbors(adata2)
sc.tl.umap(adata2)

adata = adata1.concatenate(adata2)
``` 

### 跟之前一样，你可以筛选掉那些celltypist自动注释的只有1个或者少数几个的偶然假阳性细胞
```
cell_counts = adata.obs['predicted_labels'].value_counts()
threshold = 0.005 * adata.obs.shape[0]
selected_types = cell_counts[cell_counts > threshold].index
adata = adata[adata.obs['predicted_labels'].isin(selected_types)].copy()
```

### celltypist然后concatenate后，adata仍然需要再跑一次标准流程
``` python
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata2, batch_key = 'Dataset', subset = True)  
sc.pp.scale(adata, max_value = 10)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

### 然后cellhint整合，去批次,可视化
```
cellhint.integrate(adata, 'Dataset', 'predicted_labels')

sc.tl.umap(adata)

sc.pl.umap(adata, color='predicted_labels', save='cellhint_test.png')
```

####这样你就能看到漂亮的umap分群了，adata也就是我们想要的大型单细胞图谱，我这么整合会有22万个细胞，还是很恐怖的

## 恭喜你，你可以利用自己创造的图谱进行下游的分析
















