import os
import subprocess


def process_fasta_files(input_folder, output_folder1, output_folder):
    # 确保输出文件夹存在
    os.makedirs(output_folder, exist_ok=True)

    # 遍历输入文件夹中的所有.fasta文件
    for fasta_file in os.listdir(input_folder):
        if fasta_file.endswith('.fasta'):
            base_name = os.path.splitext(fasta_file)[0]
            input_path = os.path.join(input_folder, fasta_file)
            output_path1 = os.path.join(output_folder1, f"{base_name}.fasta")
            output_path2 = os.path.join(output_folder, f"{base_name}.hmm")

            # 调用kalign命令
            subprocess.run(["kalign", "-i", input_path, "-o", output_path1], check=True)

            # 调用hmmbuild命令
            subprocess.run(["hmmbuild", output_path2, output_path1], check=True)


def combine_hmm_files(output_folder, output_file):
    # 如果输出文件已存在，则删除它
    if os.path.exists(output_file):
        os.remove(output_file)

        # 将所有.hmm文件的内容合并到输出文件中
    with open(output_file, 'wb') as outfile:
        for filename in os.listdir(output_folder):
            if filename.endswith('.hmm'):
                filepath = os.path.join(output_folder, filename)
                with open(filepath, 'rb') as infile:
                    outfile.write(infile.read())

                # 定义输入文件夹和输出文件夹


input_folder = "/home/taoyechen/NCBI_Complete_process/clu_20"
output_folder1 = "/home/taoyechen/NCBI_Complete_process/kalign_20"
output_folder = "/home/taoyechen/NCBI_Complete_process/hmm_20"
output_file = "/home/taoyechen/NCBI_Complete_process/hmm_20_combine.hmm"

# 调用函数
process_fasta_files(input_folder, output_folder1, output_folder)
combine_hmm_files(output_folder, output_file)

# 注意：这里假设所有命令都成功执行。
# 在实际使用中，你可能需要添加错误处理来捕获并处理任何失败的情况。