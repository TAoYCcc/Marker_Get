import os
import gzip
import shutil

# 源目录和目标目录
SOURCE_DIR = "/home/taoyechen/NCBI_Complete/NHP"
TARGET_DIR = "/home/taoyechen/NCBI_Complete/NHP_fna"

# 确保目标目录存在
if not os.path.exists(TARGET_DIR):
    print("目标目录不存在，正在创建...")
    os.makedirs(TARGET_DIR)

# 遍历源目录中的.gz文件
for file in os.listdir(SOURCE_DIR):
    if file.endswith('.gz'):
        # 提取文件名（不带路径和.gz扩展名）
        filename, ext = os.path.splitext(file)
        # 构建源文件的完整路径和目标文件的完整路径
        source_path = os.path.join(SOURCE_DIR, file)
        target_path = os.path.join(TARGET_DIR, filename)

        # 使用gzip模块解压文件到目标文件
        try:
            with gzip.open(source_path, 'rb') as f_in:
                with open(target_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"解压成功：{source_path} -> {target_path}")
        except Exception as e:
            print(f"解压失败：{source_path}，错误：{e}")

print("所有.gz文件已解压到", TARGET_DIR)