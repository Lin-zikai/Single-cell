# 废话不多说

## （假设model中是ENSEMBL格式，adata1中是SYMBOL）
``` python
import mygene  #没有就conda装一下
mg = mygene.MyGeneInfo()

gene_symbols = adata1.var_names.tolist()

query_results = mg.querymany(gene_symbols, scopes='symbol', fields='ensembl.gene', species='human')

symbol_to_ensembl = {}
for result in query_results:
    if 'ensembl' in result:
        ensembl_info = result['ensembl']
        if isinstance(ensembl_info, list):
            ensembl_id = ensembl_info[0].get('gene')
        else:
            ensembl_id = ensembl_info.get('gene')
        if ensembl_id:
            symbol_to_ensembl[result['query']] = ensembl_id
adata1.var_names = [symbol_to_ensembl.get(gene, gene) for gene in adata1.var_names]
```
### （假设model中是SYMBOL格式，adata1中是ENSEMBL）
``` python
import mygene
mg = mygene.MyGeneInfo()

# 假设 adata1.var_names 包含了 ENSEMBL IDs
ensembl_ids = adata1.var_names.tolist()

# 使用 mygene.querymany 查询基因符号
query_results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')

# 创建从 ENSEMBL ID 到基因符号的映射
ensembl_to_symbol = {}
for result in query_results:
    if 'symbol' in result:
        symbol = result['symbol']
        ensembl_to_symbol[result['query']] = symbol

# 更新 adata1.var_names，如果没有找到对应的基因符号则保留 ENSEMBL ID
adata1.var_names = [ensembl_to_symbol.get(ensembl_id, ensembl_id) for ensembl_id in adata1.var_names]
```
