import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import igraph as ig
import leidenalg

# 读取CSV文件
df = pd.read_csv('/home/taoyechen/parameter/allcosine_similarity_matrix.csv', index_col=0)

# 获取hmmid列表
hmmids = df.index.tolist()

# 将对角线上的值设置为1
np.fill_diagonal(df.values, 1)

# 创建CSR格式的稀疏矩阵
sparse_matrix = csr_matrix(df.values)

# 使用igraph创建图
sources, targets = sparse_matrix.nonzero()
weights = sparse_matrix[sources, targets].tolist()[0]
g = ig.Graph(directed=False)
g.add_vertices(hmmids)
g.add_edges(list(zip(sources.tolist(), targets.tolist())))
g.es['weight'] = weights

# 设置分辨率参数
resolution_parameter = 0.25  # 可以根据需要调整这个值

# 应用Leiden算法进行社区检测，包含分辨率参数
partition = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, resolution_parameter=resolution_parameter, weights='weight')

# 创建一个DataFrame来存储hmmid和它们的社区编号
community_df = pd.DataFrame({'hmmid': hmmids, 'community': partition.membership})

# 将DataFrame保存为CSV文件
community_df.to_csv('/home/taoyechen/parameter/leiden_clustering_results.csv', index=False)

# 打印操作成功的消息
print("聚类结果已保存到 'leiden_clustering_results.csv' 文件中。")
