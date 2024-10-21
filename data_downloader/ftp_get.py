import pandas as pd

# 读取 CSV 文件
df = pd.read_csv(r"C:\Users\TyCcc\Desktop\NHP.csv")  # 替换为你的 CSV 文件路径


# 定义生成 FTP 路径的函数
def generate_ftp_path(accession, name):
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    # 分割 Accession 号并构建路径
    parts = accession.split('_')
    if parts[0] == 'GCA':
        return f"{base_url}GCA/{parts[1][:3]}/{parts[1][3:6]}/{parts[1][6:9]}/{accession}_{name}"
    elif parts[0] == 'GCF':
        return f"{base_url}GCF/{parts[1][:3]}/{parts[1][3:6]}/{parts[1][6:9]}/{accession}_{name}"
    else:
        return "Invalid Accession format"

    # 应用函数到每一行


df['ftp_path'] = df.apply(lambda row: generate_ftp_path(row['Assembly Accession'], row['Assembly Name']), axis=1)

# 将修改后的 DataFrame 保存回 CSV 文件
df.to_csv(r"C:\Users\TyCcc\Desktop\NHP_ftp.csv", index=False)  # 替换为你的 CSV 文件路径
print("已保存到 CSV 文件中。")

