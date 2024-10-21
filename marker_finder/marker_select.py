import csv


def filter_csv(input_file_path, shannon_threshold, spm_threshold):
    """
    筛选CSV文件中满足特定条件的行，并将结果输出到新的CSV文件中。
    """
    with open(input_file_path, 'r', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        headers = next(csv_reader)
        selected_rows = [
            row for row in csv_reader
            if float(row[5]) > shannon_threshold and
               (float(row[3]) > spm_threshold or float(row[4]) > spm_threshold)
        ]
    return selected_rows, headers


def split_specific_clusters(filtered_rows, headers, pathogen_output_path, non_pathogen_output_path):
    """
    根据pathogen和non-pathogen的SPM值将行分配到两个不同的列表中，
    并将这些列表写入两个新的CSV文件。
    """
    pathogen_specific = []
    non_pathogen_specific = []

    for row in filtered_rows:
        if float(row[4]) > float(row[3]):
            pathogen_specific.append(row)
        else:
            non_pathogen_specific.append(row)

    with open(pathogen_output_path, 'w', newline='', encoding='utf-8') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(headers)
        csv_writer.writerows(pathogen_specific)

    with open(non_pathogen_output_path, 'w', newline='', encoding='utf-8') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(headers)
        csv_writer.writerows(non_pathogen_specific)

    print("特异蛋白簇信息已分别输出到两个文件中。")


# 初始化阈值
shannon_threshold = 0.1
spm_threshold = 0.6


# 调用filter_csv函数并获取筛选后的行和标题
filtered_rows, headers = filter_csv(
    r"C:\Users\TyCcc\Desktop\balance_tspex.csv",
    shannon_threshold,
    spm_threshold
)

# 调用split_specific_clusters函数并将筛选后的行和标题写入新的CSV文件
split_specific_clusters(
    filtered_rows,
    headers,
    'HP_marker.csv',
    'NHP_marker.csv'
)
