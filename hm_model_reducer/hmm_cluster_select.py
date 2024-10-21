import pandas as pd

# 读取CSV文件
df = pd.read_csv(r"C:\Users\TyCcc\Desktop\leiden_clustering_results.csv")

# 按community分组，并找到每个分组中shannon_specificity的最大值
max_shannon = df.groupby('community')['shannon_specificity'].transform(max)

# 创建一个布尔掩码，标记每个分组中shannon_specificity等于最大值的行
mask = df['shannon_specificity'] == max_shannon

# 在每个community分组中，如果有多个最大值，则随机选择一行
# 我们可以通过对满足条件的行进行分组并随机选择一行来实现
# 使用cumcount为每个分组中的行分配一个唯一的序号，并基于序号进行筛选
# 这里，我们假设如果有多行具有相同的最大值，我们只需要其中一行
grouped = df[mask].groupby('community', as_index=False, group_keys=False)
random_rows = []
for name, group in grouped:
    if len(group) > 1:
        # 如果有多行具有相同的最大值，随机选择一行
        random_row = group.sample(n=1)
        random_rows.append(random_row)
    else:
        # 如果只有一行具有最大值，直接添加
        random_rows.append(group)

    # 将列表中的DataFrame合并成一个DataFrame
result_df = pd.concat(random_rows, ignore_index=True)

# 将结果保存到新的CSV文件中
result_df.to_csv(r"C:\Users\TyCcc\Desktop\balance_hmm_dere.csv", index=False)

print("筛选后的数据已保存到filtered_data.csv文件中。")