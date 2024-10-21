import os
import time
import urllib
from multiprocessing import Pool, Manager

import pandas as pd
import wget

from config import download_config

config = Config()


def ftp_download(args):
    # 将需要的下载的url复制下来，左边的URL为变量，可根据实际情况填写其他的
    # args.append([real_name, fna_file_url, output_file_name, successful_list])
    name = args[0]
    fna_file_url = args[1]
    output_file_name = args[2]
    successful_list = args[3]
    e = ''
    # 检查文件是否已经下载且完整
    if os.path.exists(output_file_name) and os.path.getsize(output_file_name) > 0:
        print(f'File {output_file_name} already exists and is not empty.')
        successful_list.append(1)
        return

    start = time.time()
    for i in range(config.retry_times):
        print(f'Downloading {fna_file_url}, attempt {i + 1}/{config.retry_times}')
        try:
            # 下载文件并指定输出路径
            do = wget.download(fna_file_url, out=output_file_name)

        except urllib.error.URLError as error:
            e = error
            print(f'Error: {e}')
            now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            print(f'{now} Download of {name} failed, retrying...')
            # 移除等待时间，立即重试
            continue

        # 检查文件是否下载成功且完整
        if do and os.path.exists(output_file_name) and os.path.getsize(output_file_name) > 0:
            now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            print(f'{now} Download of {fna_file_url} successful.')
            print(f'Cost: {time.time() - start}s')
            successful_list.append(1)
            break
        else:
            print(f'Download of {fna_file_url} failed, retrying...')
            # 移除等待时间，立即重试

    # 如果文件不存在或下载失败，则记录到失败日志
    if not os.path.exists(output_file_name) or os.path.getsize(output_file_name) == 0:
        with open(config.fail_log_path, 'a') as f:
            now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            if e:
                f.write(f'{now},{fna_file_url},{e}\n')
            else:
                f.write(f'{now},{fna_file_url}\n')


def multi_download(args):
    with Pool(config.cpu_worker_num) as p:
        p.map(ftp_download, args)


def main():
    successful_list = Manager().list()
    df = pd.read_csv(config.input_path, dtype={21: str})
    ftp_list = df[config.csv_col].to_list()
    print('download list numbers:', len(ftp_list))
    args = []
    for ftp in ftp_list:
        real_name = ftp.split('/')[-1]
        fna_file_url = ftp + '/' + real_name + config.file_type1
        faa_file_url = ftp + '/' + real_name + config.file_type2
        output_file_name1 = config.output_path + real_name + config.file_type1
        output_file_name2 = config.output_path + real_name + config.file_type2
        args.append([real_name, fna_file_url, output_file_name1, successful_list])
    multi_download(args)

    # ftp_download(args[11111])


if __name__ == '__main__':
    main()
