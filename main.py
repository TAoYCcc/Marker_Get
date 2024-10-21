# main.py
import configparser
from data_downloader import download_and_extract_data
from preprocessor import preprocess_data
from hmm_model_reducer import reduce_models
from marker_finder import find_markers


def main():
    # 读取配置文件
    config = configparser.ConfigParser()
    config.read('config.ini')

    # 第一步：数据下载及解压整理
    download_and_extract_data(config['FTP']['url'], config['OUTPUT']['dir'])

    # 第二步：数据预处理
    processed_data = preprocess_data(config['OUTPUT']['dir'])

    # 第三步：模型去冗
    reduced_models = reduce_models(processed_data)

    # 第四步：特异性指标获取marker
    markers = find_markers(reduced_models, config['THRESHOLDS']['specificity'])

    # 输出结果或进行其他处理
    # ...


if __name__ == "__main__":
    main()