import gzip
import glob
import os
from multiprocessing import Pool


def process_file(file_path, output_file):
    with gzip.open(file_path, 'rb') as f:
        data = f.read()
        with open(output_file, 'ab') as out:
            out.write(data)


def merge_gz_files(input_dir, output_file):
    # 获取输入目录下所有的 .gz 文件
    files = glob.glob(os.path.join(input_dir, '*.gz'))

    # 使用多进程池来加速处理
    with Pool() as pool:
        # 对每个文件应用 process_file 函数
        for file_path in files:
            pool.apply_async(process_file, args=(file_path, output_file))

        # 关闭池并等待所有进程完成
        pool.close()
        pool.join()


# 指定输入目录和输出文件路径
input_directory = r'/home/taoyechen/NCBI_Complete/Protein'  # 替换为你的输入目录路径
output_file_path = r'/home/taoyechen/NCBI_Complete/ALL_protein.fasta'  # 替换为你的输出文件路径

# 调用函数合并文件
merge_gz_files(input_directory, output_file_path)