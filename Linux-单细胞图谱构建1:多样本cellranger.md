# 针对多个样本的cellranger，并整合(假设你所有样本的fastq都在一个文件夹）


## 
``` shell
cd /path/to/foler

for file in *fastq.gz; do
  sample_name=$(echo "$file" | awk -F '_S[0-9]+_' '{print $1}')
  if [ ! -d "$sample_name" ]; then
    mkdir "$sample_name"
  fi
  mv "$file" "$sample_name"/
done

echo "All files have been organized into their respective directories."

```

### 这样每个样本都会有一个文件夹


## 循环对每个样本运行cellranger
``` shell
export PATH=$PATH:/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/cellranger-7.2.0

TRANSCRIPTOME="/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/refdata-gex-GRCh38-2020-A"

FASTQ_ROOT="/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/HTAMSK"

for sample_dir in $FASTQ_ROOT/*; do
  if [ -d "$sample_dir" ]; then
    
    sample_name=$(basename "$sample_dir")
    
    echo "Processing $sample_name"
    
    
    output_id="${sample_name}_cellranger"
    
    cellranger count --id="$output_id" \
                     --transcriptome="$TRANSCRIPTOME" \
                     --fastqs="$sample_dir" \
                     --sample="${sample_name}" \
                    
  fi
done
```

### TRANSCRIPTOME填上参考基因组路径，FASTQ_ROOT填所有样本整合成一个个文件夹的路径

## cellranger aggr整合多个样本的cellranger结果

### 准备libs文件（样本ID,路径） （参考：https://www.jianshu.com/p/ca726a8979d7，他有一个地方有问题，就是library_id改成sample_id）
``` shell
BASE_DIR="/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/HTAMSK"

for DIR in ${BASE_DIR}/*_cellranger; do
    if [ -d "${DIR}" ]; then
        SAMPLE_ID=$(basename ${DIR})
        H5_PATH="${DIR}/outs/molecule_info.h5"
        if [ -f "${H5_PATH}" ]; then
            echo "${SAMPLE_ID},${H5_PATH}" >> libs.csv
        else
            echo "Warning: File not found - ${H5_PATH}"
        fi
    fi
done
```

###最后运行
``` shell
cellranger aggr --id=HTAMSK_multi_sample --csv=libs.csv --normalize=mapped
```

### --id是输出文件夹的名字，--csv是libs文件，最后输出一个文件夹，里面的outs/counts就是整合好的单细胞了

##恭喜你完成了多个样本的整合，下面我们第二个是分享，批量分群+去批次构建大型单细胞图谱











