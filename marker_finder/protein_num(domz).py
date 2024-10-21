import os
import re
import csv
from collections import defaultdict

# 假设这是包含hmmid的txt文件路径
hmmid_file_path = '/home/taoyechen/Balance_process/hmmid.txt'
# 假设这是包含比对结果文件夹的父目录路径
alignment_folder_base_path = '/home/taoyechen/Balance_process/nhp_result'
# CSV文件输出路径
csv_output_path = '/home/taoyechen/Balance_process/nhp_result_domZ.csv'

# 读取hmmid文件并获取所有的hmmid
with open(hmmid_file_path, 'r') as hmmid_file:
    hmmids = [line.strip() for line in hmmid_file.readlines()]

# 初始化一个空列表来存储所有需要处理的文件路径
filenames = []

# 遍历所有的hmmid，并收集对应的比对结果文件路径
for hmmid in hmmids:
    # 构造hmmid对应的文件夹路径
    hmmid_folder_path = os.path.join(alignment_folder_base_path, hmmid)

    # 检查文件夹是否存在
    if os.path.exists(hmmid_folder_path) and os.path.isdir(hmmid_folder_path):
        # 遍历文件夹中的txt文件（比对结果文件）
        for alignment_file in os.listdir(hmmid_folder_path):
            if alignment_file.endswith('.txt'):
                filenames.append(os.path.join(hmmid_folder_path, alignment_file))

file_count = len(filenames)
print(f"文件路径个数: {file_count}")

# 检查filenames列表的长度，确保不会超出索引范围
if file_count >= 10:
    # 使用切片查看前十条文件路径
    first_ten_files = filenames[:10]
    for file_path in first_ten_files:
        print(file_path)
else:
    # 如果filenames列表中的文件路径少于10条，则打印所有文件路径
    print("Files found:")
    for file_path in filenames:
        print(file_path)

domz_values = {}  # 初始化字典以存储 domZ 值

for file_path in filenames:
    # 解析文件名以获取 hmmid 和 sequence_id
    file_name = os.path.basename(file_path)
    parts = file_name.split('_')

    # 根据您的规则提取 hmmid 和 sequence_id
    if len(parts) > 5:  # 确保文件名有足够的部分
        hmmid = '_'.join(parts[0:2])  # 第二个分隔符前的所有字符，即 parts[1] 和 parts[2]
        sequence_id = '_'.join(parts[3:5])  # 第三到第五个分隔符间的所有字符，即 parts[3], parts[4], 和 parts[5]
        # print(f"hmmid: {hmmid}, sequence_id: {sequence_id}")
    else:
        pass
        # print("Filename does not match the expected format.")

    # 读取文件并提取 domZ 值
    with open(file_path, 'r') as file:
        data_lists = [line.strip() for line in file.readlines()]
        data = data_lists[-5]
        new_data1 = data.split(':')[1]
        domZ_value = new_data1.split('[')[0].strip()
        # 将 domZ 值存储在字典中，键是 (hmmid, sequence_id) 的元组
        domz_values[(hmmid, sequence_id)] = domZ_value
    # text = "hmmid: {}, sequence_id: {}, domZ_value: {}".format(hmmid, sequence_id, domZ_value)

# 提取所有独特的hmmid和sequence_id
unique_hmmids = sorted(list(set([key[0] for key in domz_values.keys()])))
unique_sequence_ids = sorted(list(set([key[1] for key in domz_values.keys()])))

# 在写入CSV之前，打印domz_values以确认数据存在
# print("domz_values:", domz_values)

# 打印unique_hmmids和unique_sequence_ids以确认它们不为空
# print("unique_hmmids:", unique_hmmids)
# print("unique_sequence_ids:", unique_sequence_ids)
# 打开CSV文件进行写入
with open(csv_output_path, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)

    # 写入CSV的头部，即所有的sequence_id
    header = ['hmmid'] + unique_sequence_ids
    writer.writerow(header)

    # 遍历所有的hmmid，并为每个hmmid写入一行数据
    for hmmid in unique_hmmids:
        row = [hmmid]  # 第一列是hmmid
        for sequence_id in unique_sequence_ids:
            # 查找对应的domZ值，如果不存在则填充None或空字符串
            domZ_value = domz_values.get((hmmid, sequence_id), '')
            row.append(domZ_value)
        writer.writerow(row)  # 写入CSV文件

print(f"CSV文件已保存到 {csv_output_path}")
