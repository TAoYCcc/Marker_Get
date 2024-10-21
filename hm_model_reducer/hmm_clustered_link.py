import os

# 设置目录和文件路径
hmmid_file = "/home/taoyechen/Balance_process/hmmid.txt"
hmm_models_dir = "/home/taoyechen/Balance_process/hmm_20"
links_dir = "/home/taoyechen/Balance_process/hmm_clustered"

# 确保链接目录存在
os.makedirs(links_dir, exist_ok=True)

# 读取hmmid.txt中的每一行
with open(hmmid_file, 'r') as file:
    for hmmid in file:
        hmmid = hmmid.strip()  # 去除行尾的换行符
        # 构建hmm模型文件的完整路径
        hmm_file = os.path.join(hmm_models_dir, f"{hmmid}.hmm")
        # 构建链接的完整路径
        link_file = os.path.join(links_dir, f"{hmmid}.hmm")

        # 检查hmm模型文件是否存在
        if os.path.isfile(hmm_file):
            # 如果链接文件已存在，则删除它（如果需要的话）
            if os.path.exists(link_file):
                os.remove(link_file)
                # 创建到hmm模型文件的符号链接
            os.symlink(hmm_file, link_file)
            print(f"Created link: {link_file} -> {hmm_file}")
        else:
            print(f"No hmm file found for {hmmid}")

print("All links have been created or checked.")