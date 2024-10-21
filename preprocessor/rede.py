import subprocess
import os

# 设置输入目录和MMseqs2数据库路径
input_dir = "/home/taoyechen/NCBI_Complete_process"
dr = os.path.join(input_dir, "dere/dr_db")
drclu = os.path.join(input_dir, "dere/drc_db")
tmp_dir = os.path.join(input_dir, "tmp")
tsv_file = os.path.join(input_dir, "dere/drclu.tsv")
dbrep = os.path.join(input_dir, "dere/DB_clu_rep")
fasta_output = os.path.join(input_dir, "clusters_rep.fasta")
all_clusters_fasta = os.path.join(input_dir, "all_clusters.fasta")

# 创建MMseqs2数据库
subprocess.run(["mmseqs", "createdb", os.path.join(input_dir, "ALL_protein.fasta"), dr], check=True)

# 使用MMseqs2的linclust对蛋白序列进行去冗余
subprocess.run([
    "mmseqs", "linclust", dr, drclu, tmp_dir,
    "--threads", "150", "--kmer-per-seq", "80", "-c", "1.0",
    "--cluster-mode", "2", "--cov-mode", "1", "--min-seq-id", "0.95"
], check=True)

# 创建tsv文件
subprocess.run(["mmseqs", "createtsv", dr, dr, drclu, tsv_file], check=True)

# 提取每个cluster的代表序列
subprocess.run(["mmseqs", "createsubdb", drclu, dr, dbrep], check=True)
subprocess.run(["mmseqs", "convert2fasta", dbrep, fasta_output], check=True)

# Extract all cluster sequences
subprocess.run(["mmseqs", "createseqfiledb", dr, drclu, "DB_clu_seq"], check=True)
subprocess.run(["mmseqs", "result2flat", dr, dr, "DB_clu_seq", all_clusters_fasta], check=True)

# 注意：这里假设所有MMseqs2命令都成功执行。
# 在实际使用中，你可能需要添加错误处理来捕获并处理任何失败的情况。