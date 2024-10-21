import os


# 读取特异蛋白簇ID的函数
def read_ids_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]


def combine_models(ids_list, models_directory, combined_model_path):
    # 打开合并文件，准备写入数据
    with open(combined_model_path, 'ab') as combined_file:
        for model_id in ids_list:
            model_name = model_id + '.hmm'  # 假设模型文件的扩展名为.hmm
            model_path = os.path.join(models_directory, model_name)

            # 检查模型文件是否存在
            if os.path.isfile(model_path):
                # 读取模型文件内容
                with open(model_path, 'rb') as model_file:
                    model_content = model_file.read()
                # 将模型内容追加到合并文件中
                combined_file.write(model_content)
                print(f"模型 {model_name} 的内容已追加到 {combined_model_path} 中")
            else:
                print(f"模型 {model_name} 不存在于 {models_directory}")


# 特异蛋白簇ID列表的文件路径
pathogen_ids_file = '/home/taoyechen/classification_pre/New_feature_get/pm.txt'  # 替换为pathogen特异蛋白簇ID列表的实际文件路径
non_pathogen_ids_file = '/home/taoyechen/classification_pre/New_feature_get/nonm.txt'  # 替换为non_pathogen特异蛋白簇ID列表的实际文件路径

# 包含HMM模型的文件夹路径
models_directory = '/home/taoyechen/classification_pre/hmm_20'  # 替换为HMM模型文件夹的实际路径

# 目标目录路径
nm_model_path = '/home/taoyechen/classification_pre/New_feature_get/nonm.hmm'  # 替换为pathogen特异模型的目标目录路径
pm_model_path = '/home/taoyechen/classification_pre/New_feature_get/pm.hmm'  # 替换为pathogen特异模型的目标目录路径

# 读取ID并创建符号链接
pathogen_ids = read_ids_from_file(pathogen_ids_file)
non_pathogen_ids = read_ids_from_file(non_pathogen_ids_file)
combine_models(pathogen_ids, models_directory, nm_model_path)
combine_models(non_pathogen_ids, models_directory, pm_model_path)

