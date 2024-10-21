import pandas as pd

# 读取文本文件内容，跳过注释行
data = []
with open(r'C:\Users\TyCcc\Desktop\test1.txt', 'r') as file:
    for line in file:
        if not line.startswith('#') and line.strip():
            data.append(line.split())

# 定义列名
columns = [
    'seq_id', 'accession', 'hmm_id', 'query_accession',
    'full_E-value', 'full_score', 'full_bias', 'best_1_domain_E-value',
    'best_1_domain_score', 'best_1_domain_bias', 'exp', 'reg', 'clu',
    'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target'
]

# 创建DataFrame
df = pd.DataFrame(data, columns=columns)

# 将DataFrame保存到CSV文件
df.to_csv(r'C:\Users\TyCcc\Desktop\output.csv', index=False)