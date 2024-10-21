import os
import pandas as pd


def process_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        # 跳过标题行
        header = next(file).strip().split('\t')
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 7 and parts[1] == 'rRNA':
                try:
                    placement = int(parts[6])
                    placement_div_3 = placement / 3
                except ValueError:
                    # 如果Placements不是数字（例如'na'），则跳过该行或设置为NaN
                    placement_div_3 = float('nan')  # 或者你可以选择跳过这行数据：continue
                data.append([parts[3], placement_div_3])  # Assembly-unit accession 和 Placements数除以3的结果
    return data, header


def process_folder(folder_path):
    all_data = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(folder_path, filename)
            file_data, header = process_file(file_path)
            all_data.extend(file_data)
    return all_data, header


def save_to_excel(data, header, output_path):
    # 创建一个DataFrame，但我们需要手动设置列名，因为header中的列名可能不完全匹配我们的数据
    df = pd.DataFrame(data, columns=['Assembly-unit accession', 'Placements / 3'])
    # 如果需要，可以将原始的header信息保存到Excel的一个额外sheet或作为注释
    # 这里我们仅保存处理后的数据
    df.to_excel(output_path, index=False)


if __name__ == "__main__":
    folder_path = "/home/taoyechen/NCBI_oneS_balance/HP_rna_txt/" # 替换为你的文件夹路径
    output_excel_path = '/home/taoyechen/NCBI_oneS_balance/HP_data.xlsx'  # 输出Excel文件的路径
    data, _ = process_folder(folder_path)
    # 由于我们重新定义了列名，所以不需要原始的header信息来创建DataFrame
    # 但如果你需要它，可以在process_folder函数中保留并传递它
    save_to_excel(data, None, output_excel_path)
    print(f"处理后的数据已保存到 {output_excel_path}")
