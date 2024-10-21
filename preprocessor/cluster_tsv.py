import subprocess


def create_tsv_file(input_file, db_one, db_two, tsv_file):
    # 创建MMseqs2数据库
    subprocess.run(["mmseqs", "createdb", input_file, db_one], check=True)

    # 使用MMseqs2的cluster对去冗余后的蛋白序列进行去聚类
    subprocess.run([
        "mmseqs", "cluster", db_one, db_two, "tmp",
        "-s", "4", "-e", "1e-5", "-c", "0.8", "--cov-mode", "0", "--cluster-mode", "0",
        "--max-seqs", "5000", "--min-seq-id", "0.5", "--cluster-reassign", "1"
    ], check=True)

    # 创建tsv文件
    subprocess.run(["mmseqs", "createtsv", db_one, db_one, db_two, tsv_file], check=True)


# 指定输入文件、输出文件以及MMseqs2数据库的路径
input_file = "/home/taoyechen/NCBI_Complete_process/clusters_rep.fasta"
output_file = "/home/taoyechen/NCBI_Complete_process/cluster/clustered.tsv"
one = "/home/taoyechen/NCBI_Complete_process/cluster/rede_db"
two = "/home/taoyechen/NCBI_Complete_process/cluster/clu_db"

# 调用函数
create_tsv_file(input_file, one, two, output_file)

# 注意：这里假设所有MMseqs2命令都成功执行。
# 在实际使用中，你可能需要添加错误处理来捕获并处理任何失败的情况。