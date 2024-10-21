import os
import csv
from multiprocessing import Pool

# 指定目录路径
directory = '/home/taoyechen/classification_pre/hmm_align/path/result/'

# 设置文件大小阈值（以字节为单位）
size_threshold = 2 * 1024   # 2kb

# 定义一个函数来计算超过阈值的文件数量
def count_files(folder):
    count = 0
    for filename in os.listdir(folder):
        filepath = os.path.join(folder, filename)
        if os.path.isfile(filepath):
            if os.path.getsize(filepath) > size_threshold:
                count += 1
    return count

# 使用Pool来创建指定数量的进程
def process_folders(directory, processes=20):
    folders = [os.path.join(directory, folder) for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder))]
    with Pool(processes) as pool:
        results = pool.map(count_files, folders)
    return results

# 运行函数并将结果写入CSV文件
if __name__ == '__main__':
    results = process_folders(directory)
    csv_file_path = '/home/taoyechen/classification_pre/hmm_align/path_counts.csv'
    with open(csv_file_path, 'w', newline='') as csvfile:
        fieldnames = ['Folder', 'Files_Above_Threshold']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for folder, count in zip(os.listdir(directory), results):
            writer.writerow({'Folder': folder, 'Files_Above_Threshold': count})
    print(f'计数完成，结果已输出到folder_file_counts.csv文件中。')