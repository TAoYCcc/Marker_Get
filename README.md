# 数据下载及整理

## Ftp_get

```python
import pandas as pd

# 读取 Excel 文件
df = pd.read_excel(r"C:\Users\TyCcc\Desktop\3\ncbi_5_8.xlsx")  # 替换为你的 Excel 文件路径

# 在将'Assembly Name'传递给函数之前替换空格为下划线  
df['Assembly Name'] = df['Assembly Name'].str.replace(' ', '_') 

# 定义生成 FTP 路径的函数
def generate_ftp_path(accession, name):
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    # 分割 Accession 号并构建路径
    parts = accession.split('_')
    if parts[0] == 'GCA':
        return f"{base_url}GCA/{parts[1][:3]}/{parts[1][3:6]}/{parts[1][6:9]}/{accession}_{name}"
    elif parts[0] == 'GCF':
        return f"{base_url}GCF/{parts[1][:3]}/{parts[1][3:6]}/{parts[1][6:9]}/{accession}_{name}"
    else:
        return "Invalid Accession format"

# 应用函数到每一行
df['ftp_path'] = df.apply(lambda row: generate_ftp_path(row['Assembly Accession'], row['Assembly Name']), axis=1)

# 将修改后的 DataFrame 保存回 Excel 文件
df.to_excel(r"C:\Users\TyCcc\Desktop\3\ncbi_58_ftp.xlsx", index=False)  # 替换为你的 Excel 文件路径
print("已保存到 Excel 文件中。")

```

## NCBI_download

### Download_main

```python
import os
import time
import urllib
from multiprocessing import Pool, Manager

import pandas as pd
import wget

from download_config import Config

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

```

### Download_config

```python
class Config:
    def __init__(self):
        self.input_path = r'/home/taoyechen/NCBI_Complete/NHP.csv'
        self.output_path = r'/home/taoyechen/NCBI_Complete/NHP/'
        self.cpu_worker_num = 180
        # file suffix of final url
        self.file_type1 = '_genomic.fna.gz'
        self.file_type2 = '_protein.faa.gz'
        self.csv_col = 'ftp_path'

        # times of re-downloads
        self.retry_times = 10
        self.fail_log_path = r'/home/taoyechen/NCBI_Complete/NHP_log.txt'
```

## Merge_gz2allfaa

```python
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

```



## Unzip

### Unzip4zip

```shell
#!/bin/bash


ZIP_FILES=$(ls *.zip) #获取当前目录下所有.zip结尾的文件


for zip_file in $ZIP_FILES; do
        fname=$(basename $zip_file .zip)
        echo "$fname"
        ZIP_TO="/home/taoyechen/predata/pathogen/raw_data/$fname"
        mkdir /home/taoyechen/predata/pathogen/raw_data/$fname
        # 开始解压
        #[注:  x  参数解压；
        #     -o 参数解压文件位置；     ]
        7z x $zip_file -r -o$ZIP_TO;

done
```

### Unzip4gz

```shell
#!/bin/bash  
  
# 源目录，即包含.gz文件的目录  
SOURCE_DIR="/home/taoyechen/NCBI_oneS_balance/NHP_rna" 
# 目标目录，即你想要解压到的目录  
TARGET_DIR="/home/taoyechen/NCBI_oneS_balance/NHP_rna_tx"  
  
# 确保目标目录存在  
if [ ! -d "$TARGET_DIR" ]; then  
    echo "目标目录不存在，正在创建..."  
    mkdir -p "$TARGET_DIR"  
fi  
  
# 遍历源目录中的.gz文件  
for file in "$SOURCE_DIR"/*.gz; do  
    # 提取文件名（不带路径和.gz扩展名）  
    filename=$(basename "$file" .gz)  
    # 构建目标文件的路径  
    target_file="$TARGET_DIR/$filename"  
    # 使用gunzip -c解压文件到目标文件  
    gunzip -c "$file" > "$target_file"  
    # 检查是否解压成功  
    if [ $? -eq 0 ]; then  
        echo "解压成功：$file -> $target_file"  
    else  
        echo "解压失败：$file"  
    fi  
done  
  
echo "所有.gz文件已解压到$TARGET_DIR"
```



## Link_faa/fna

```shell
#!/bin/bash

# 指定目录，替换成想要查找文件的目录
specified_directory="/path/to/your/specified/directory"
# 目标目录，替换成想要链接文件的目录
target_directory="/home/taoyechen/classification_pre/path/protein_seq"


# 查找指定目录下的所有以.faa结尾的文件
find "$specified_directory" -type f -name "*.faa" | while read filename; do
  # 获取文件的母文件夹名称
  parent_folder="$(basename "$(dirname "$filename")")"

  # 构建新的文件名，格式为母文件夹名称.faa
  new_filename="${parent_folder}.faa"

  # 在目标目录下创建链接
  ln -s "$(realpath "$filename")" "$target_directory/$new_filename"
done
```

## Combine_2allfaa

```shell
#!/bin/bash

# 指定输入目录和输出目录
input_dir="/home/taoyechen/classification_pre/mmseqs2_ctr/protein_seq"
output_dir="/home/taoyechen/classification_pre/mmseqs2_ctr"

# 检查输出目录是否存在，如果不存在则创建
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

# 指定输出文件
output_file="$output_dir/combined_proteins.fasta"

# 检查输出文件是否存在，如果存在则删除
if [ -f "$output_file" ]; then
  rm "$output_file"
fi

# 遍历输入目录下的所有.faa文件
for input_file in "$input_dir"/*.faa; do
  # 将输入文件的内容追加到输出文件
  cat "$input_file" >> "$output_file"
done
```



# Data_process

## Predict_protein(outofdate)

```shell
#!/bin/bash

# 指定输入和输出目录
input_dir="/home/taoyechen/NCBI_Complete/HP_fna"
output_dir="/home/taoyechen/NCBI_Complete/HP_faa"

# 检查输出目录是否存在，如果不存在则创建
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

# 遍历输入目录下的所有.fna文件
for input_file in "$input_dir"/*.fna; do
  # 获取输入文件的基本名（不包括路径和扩展名）
  base_name=$(basename "$input_file" .fna)
  # 创建输出文件的路径
  output_file="$output_dir/${base_name}_prodigal.gff"
  # 运行Prodigal软件
  prodigal -i "$input_file" -o "$output_file" -a "${output_file/.gff/.faa}" -f gff
done
```

## Rede

### Protein_rede

```shell
#!/bin/bash

input_dir="/home/taoyechen/NCBI_Complete_process"

# MMseqs2数据库
dr="$input_dir/dere/dr_db" 

# 创建MMseqs2数据库
mmseqs createdb "$input_dir/ALL_protein.fasta" "$dr"

# 指定去冗余后的数据库
drclu="$input_dir/dere/drc_db" 

# 使用MMseqs2的linclust对蛋白序列进行去冗余
mmseqs linclust "$dr" "$drclu" "$input_dir/tmp" --threads 150 --kmer-per-seq 80 -c 1.0 --cluster-mode 2 --cov-mode 1 --min-seq-id 0.95

# 创建tsv文件
mmseqs createtsv "$dr" "$dr" "$drclu" "$input_dir/dere/drclu.tsv"

# 指定去冗余后的代表序列数据库
dbrep="$input_dir/dere/DB_clu_rep" 

# 提取每个 cluster 的代表序列
mmseqs createsubdb "$drclu" "$dr" "$dbrep"
mmseqs convert2fasta "$dbrep" clusters_rep.fasta

# Extract fasta
mmseqs createseqfiledb "$dr" "$drclu" DB_clu_seq
mmseqs result2flat "$dr" "$dr" DB_clu_seq all_clusters.fasta
```

### Example_code

```shell
# Place of fasta and database
input_faa=/some-place/input.faa
mmseq_db_home=/some-place

DB_name=prot90
DB_path=${mmseq_db_home}/${DB_name}

min_cov=0.8 # 至少 80% 的 coverage
min_id=0.9 # 至少 90% 的 sequence identity

# 1. fasta to DB format
mkdir ${DB_path} # 此时创建了 /some-place/prot90 这个文件夹，但是此时的 prot90 文件夹为空
tmp_folder=${DB_path}/tmp
mkdir ${tmp_folder} 
DB=${DB_path}/${DB_name} 

# mmseqs createdb proteome.fasta DB
mmseqs createdb ${input_faa} ${DB}  # 这一步会在 prot90 这个文件夹里创建 prot.xxx 的 DB 文件

# 2. Clustering
DB_clu_folder=${DB_path}/clu
mkdir ${DB_clu_folder}
DB_clu=${DB_clu_folder}/clu

# mmseqs cluster DB DB_clu tmp
mmseqs cluster ${DB} ${DB_clu} ${tmp_folder} --cov-mode 0 -c ${min_cov} --min-seq-id ${min_id} --threads 8 

# 3. Extract tsv
# mmseqs createtsv DB DB DB_clu clustered.tsv
mmseqs createtsv ${DB} ${DB} ${DB_clu} ${DB_path}/clustered.tsv # 记录有cluster 的信息
 
# 4. Extract fasta
# mmseqs createseqfiledb DB DB_clu DB_clu_seq
# mmseqs result2flat DB DB DB_clu_seq DB_clu_seq.fasta --use-fasta-header
DB_clu_seq_folder=${DB_path}/clu_seq
mkdir ${DB_clu_seq_folder}
DB_clu_seq=${DB_clu_seq_folder}/clu_seq

mmseqs createseqfiledb ${DB} ${DB_clu} ${DB_clu_seq}
mmseqs result2flat ${DB} ${DB} ${DB_clu_seq} ${DB_clu_seq_folder}/clu_seq.fasta --use-fasta-header # 记录cluster 内部的序列信息
 
# 5. Extract representative sequence 
# mmseqs createsubdb DB_clu DB DB_clu_rep
# mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta
DB_clu_rep_folder=${DB_path}/clu_rep
mkdir ${DB_clu_rep_folder}
DB_clu_rep=${DB_clu_rep_folder}/clu_rep

mmseqs createsubdb ${DB_clu} ${DB} ${DB_clu_rep}
mmseqs convert2fasta ${DB_clu_rep} ${DB_clu_rep_folder}/clu_rep.fasta # 代表序列信息
```

### Rede_all_protein(python_by_tsv)

```python
import os

# 替换为您的TSV文件路径
tsv_file = '/home/taoyechen/classification_pre/mmseqs2_ctr/dere/drclu.tsv'
# 替换为您的fasta文件路径
fasta_file = '/home/taoyechen/classification_pre/mmseqs2_ctr/combined_proteins.fasta'
# 输出文件夹
output_dir = '/home/taoyechen/classification_pre/mmseqs2_ctr/all_clusters.fasta'

# 读取fasta文件并创建一个字典，键为序列ID，值为序列
sequences = {}
with open(fasta_file, 'r') as fasta:
    sequence_id = ''
    sequence = ''
    for line in fasta:
        if line.startswith('>'):
            if sequence_id:
                sequences[sequence_id] = sequence
            sequence_id = line.strip().split()[0][1:]  # 移除'>'并取第一部分作为ID
            sequence = ''
        else:
            sequence += line.strip()
    if sequence_id:
        sequences[sequence_id] = sequence

# 读取TSV文件并将每个聚类的序列合并到一个文件中
with open(tsv_file, 'r') as tsv, open(output_dir, 'w') as output_fasta:
    for line in tsv:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            cluster_id, sequence_id = parts
            output_fasta.write(f'>{sequence_id}\n{sequences[sequence_id]}\n')
```

## Cluster

### Cluster_tsv

```shell
#!/bin/bash


create_tsv_file() {
  local input_file=$1
  local db_one=$2
  local db_two=$3
  local tsv_file=$4

  # 创建MMseqs2数据库
  mmseqs createdb "$input_file" "$db_one"

  # 使用MMseqs2的clust对去冗余后的蛋白序列进行去聚类
  mmseqs cluster "$db_one" "$db_two" tmp -s 4 -e 1e-5 -c 0.8 --cov-mode 0 --cluster-mode 0 --max-seqs 5000 --min-seq-id 0.5 --cluster-reassign 1

  # 创建tsv文件
  mmseqs createtsv "$db_one" "$db_one" "$db_two" "$tsv_file"
}


input_file="/home/taoyechen/NCBI_Complete_process/clusters_rep.fasta"
output_file="/home/taoyechen/NCBI_Complete_process/cluster/clutered.tsv"
# 指定需要处理的MMseqs2数据库
one="/home/taoyechen/NCBI_Complete_process/cluster/rede_db"
# 指定去冗余后的数据库
two="/home/taoyechen/NCBI_Complete_process/cluster/clu_db"

create_tsv_file "$input_file" "$one" "$two" "$output_file"

```

### Cluster_allclusters_by_mmseqs2

```shell
#!/bin/bash

input_dir="/home/taoyechen/NCBI_Complete_process"

# MMseqs2数据库
dr="$input_dir/cluster/dere_db" 

# 创建MMseqs2数据库
mmseqs createdb "$input_dir/clusters_rep.fasta" "$dr"

# 指定聚类后的数据库
drclu="$input_dir/dere/clu_db" 

# 使用MMseqs2的clust对去冗余后的蛋白序列进行去聚类
  mmseqs cluster "$dr" "$drclu" "$input_dir/tmp" --thread 150 -s 4 -e 1e-5 -c 0.8 --cov-mode 0 --cluster-mode 0 --max-seqs 5000 --min-seq-id 0.5 --cluster-reassign 1 
  


# 创建tsv文件
mmseqs createtsv "$dr" "$dr" "$drclu" "$input_dir/cluster/clu.tsv"

# Extract fasta
mmseqs createseqfiledb "$dr" "$drclu" DB_clu_seq
mmseqs result2flat "$dr" "$dr" DB_clu_seq all_clusters.fasta
```

### Cluster_protein_20(python_by_tsv)

```python
import os

def read_fasta(fasta_file):
    """Read a FASTA file and return a dictionary with sequence IDs as keys and sequences as values."""
    sequences = {}
    with open(fasta_file, 'r') as fasta:
        sequence_id = ''
        sequence = ''
        for line in fasta:
            if line.startswith('>'):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line.strip().split()[0][1:]  # Remove '>' and take the first part as ID
                sequence = ''
            else:
                sequence += line.strip()
        if sequence_id:
            sequences[sequence_id] = sequence
    return sequences

def write_clusters(tsv_file, sequences, output_dir, threshold):
    """Read a TSV file and write sequences of each cluster to separate files if they exceed the threshold."""
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cluster_counts = {}
    with open(tsv_file, 'r') as tsv:
        cluster_sequences = {}
        for line in tsv:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                cluster_id, sequence_id = parts
                if cluster_id not in cluster_sequences:
                    cluster_sequences[cluster_id] = []
                cluster_sequences[cluster_id].append(sequence_id)

    for cluster_id, sequence_ids in cluster_sequences.items():
        cluster_count = len(sequence_ids)
        cluster_counts[cluster_id] = cluster_count
        if cluster_count > threshold:
            with open(os.path.join(output_dir, f'cluster_{cluster_id}.fasta'), 'w') as cluster_file:
                for sequence_id in sequence_ids:
                    cluster_file.write(f'>{sequence_id}\n{sequences[sequence_id]}\n')
    return cluster_counts


tsv_file = '/home/taoyechen/classification_pre/mmseqs2_ctr/clu_ed/fin.tsv'
fasta_file = '/home/taoyechen/classification_pre/mmseqs2_ctr/all_clusters.fasta'
output_dir = '/home/taoyechen/classification_pre/mmseqs2_ctr/clu_protein_20'
threshold = 20  

# Read the fasta file and create a dictionary of sequences
sequences = read_fasta(fasta_file)

# Write cluster sequences to separate files and get the counts
cluster_counts = write_clusters(tsv_file, sequences, output_dir, threshold)
```

### Cluster_protein_20(shell_by_tsv)

```shell
#!/bin/bash

# 设置阈值
threshold=20  # 根据需要更改这个阈值

read_fasta() {
  local fasta_file=$1
  declare -A sequences
  while read -r line; do
    if [[ $line == '>'* ]]; then
      if [[ -n $sequence_id ]]; then
        sequences[$sequence_id]=$sequence
      fi
      sequence_id=${line#>} # 移除'>'并取第一部分作为ID
      sequence_id=${sequence_id%% *} # 取空格前的部分作为ID
      sequence=''
    else
      sequence+=$line
    fi
  done < "$fasta_file"
  if [[ -n $sequence_id ]]; then
    sequences[$sequence_id]=$sequence
  fi
  echo ${sequences[@]}
}

write_clusters() {
  local tsv_file=$1
  local output_dir=$2
  local sequences
  IFS=' ' read -r -a sequences <<< "$3"
  mkdir -p "$output_dir"
  declare -A cluster_counts
  while IFS=$'\t' read -r cluster_id sequence_id; do
    if [[ -n ${sequences[$sequence_id]} ]]; then
      if [[ ! -v cluster_counts[$cluster_id] ]]; then
        cluster_counts[$cluster_id]=0
      fi
      ((cluster_counts[$cluster_id]++))
      if [[ ${cluster_counts[$cluster_id]} -gt $threshold ]]; then
        echo ">$sequence_id" >> "$output_dir/cluster_${cluster_id}.fasta"
        echo "${sequences[$sequence_id]}" >> "$output_dir/cluster_${cluster_id}.fasta"
      fi
    fi
  done < "$tsv_file"
  for cluster_id in "${!cluster_counts[@]}"; do
    echo "Cluster $cluster_id has ${cluster_counts[$cluster_id]} sequences."
  done
}

# 替换为您的TSV文件路径
tsv_file='/path/to/your/tsv_file.tsv'
# 替换为您的fasta文件路径
fasta_file='/path/to/your/fasta_file.fasta'
# 输出文件夹
output_dir='/path/to/your/output_directory'

# 读取fasta文件并创建序列数组
sequences=$(read_fasta "$fasta_file")

# 将每个簇的序列写入到单独的文件中，并统计数量
write_clusters "$tsv_file" "$output_dir" "$sequences"
```

## Align_hmmbu

```shell
#!/bin/bash

process_fasta_files() {
  local input_folder=$1
  local output_folder1=$2
  local output_folder=$3

  mkdir -p $output_folder

  for fasta_file in $input_folder/*.fasta; do
    local base_name=$(basename $fasta_file .fasta)
    kalign -i $fasta_file -o $output_folder1/${base_name}.fasta
    hmmbuild $output_folder/${base_name}.hmm $output_folder1/${base_name}.fasta
  done
}

combine_hmm_files() {
  local output_folder=$1
  local output_file=$2

  if [ -f "$output_file" ]; then
    rm "$output_file"
  fi

  for file in "$output_folder"/*.hmm; do
    cat "$file" >> "$output_file"
  done
}

# 定义输入文件夹和输出文件夹
input_folder="/home/taoyechen/NCBI_Complete_process/clu_20"
output_folder1="/home/taoyechen/NCBI_Complete_process/kalign_20"
output_folder="/home/taoyechen/NCBI_Complete_process/hmm_20"
output_file="/home/taoyechen/NCBI_Complete_process/hmm_20_combine.hmm"

# 调用函数
process_fasta_files "$input_folder" "$output_folder1" "$output_folder"
combine_hmm_files "$output_folder" "$output_file"


```

## Hmm_align

### Hmm_align(allseq)

```shell
#!/bin/bash


# 定义HMM模型文件夹和序列数据库文件夹的路径
hmm_folder="/home/taoyechen/Balance_process/hmm_clustered"
sequence_db_folder="/home/taoyechen/NCBI_oneS_balance/HP_faa_296"
res_folder="/home/taoyechen/Balance_process/HP2hmmC"


# 设置并发数
concurrency=20
current_jobs=0

# 循环处理每个HMM模型文件
for hmm_file in $hmm_folder/*.hmm; do
    # 获取不包含路径和扩展名的文件名
    base_name=$(basename $hmm_file .hmm)

    # 为每个模型创建一个结果文件夹
    mkdir -p "$res_folder/$base_name"

    # 循环处理每个序列数据库文件
    for sequence_db in $sequence_db_folder/*.faa; do
        # 获取序列数据库文件的基本名
        db_name=$(basename $sequence_db .faa)

        # 运行hmmsearch并将结果输出到相应的文件夹中
        (hmmsearch --cpu 10 -E 1e-5 $hmm_file $sequence_db > "$res_folder/$base_name/${base_name}_to_${db_name}_results.txt") &

    # 更新当前作业计数
        ((current_jobs++))

        # 如果达到并发数，等待所有后台作业完成
        if ((current_jobs >= concurrency)); then
            wait
            current_jobs=0
        fi
    done
done


# 等待最后一批作业完成
wait
```

### Hmm_align(single_seq)

```shell
#!/bin/bash  
  
# HMM模型文件所在的目录  
HMM_DIR="/home/taoyechen/Balance_process/hmm_clustered"  
# Fasta文件所在的目录  
FAA_DIR="/home/taoyechen/NCBI_oneS_balance/GCA_001499615.1_Acetobacter_senegalensis_108B_protein.faa"  
# 输出结果的目录  
OUTPUT_DIR="/home/taoyechen/Balance_process/add"  
  
# 确保输出目录存在  
mkdir -p "$OUTPUT_DIR"  
  
# 遍历HMM模型目录  
for hmm_file in "$HMM_DIR"/*.hmm; do  
    # 提取HMM模型的ID（这里假设文件名是hmmid.hmm）  
    hmm_id=$(basename -- "$hmm_file" .hmm)  
  
    # 遍历Fasta文件目录  
    for faa_file in "$FAA_DIR"/*.faa; do  
        # 提取Fasta文件名（不带扩展名）  
        faa_name=$(basename -- "$faa_file" .faa)  
  
        # 构建输出文件名  
        output_file="$OUTPUT_DIR/${hmm_id}_to_${faa_name}.out"  
  
        # 使用hmmsearch进行比对，并将结果保存到输出文件  
        hmmsearch --cpu 10 -E 1e-5 "$hmm_file" "$faa_file" > "$output_file"  
  
        # 输出信息到控制台（可选）  
        echo "Processed $hmm_file to $faa_file and saved to $output_file"  
    done  
done
```

## Count_gemone

```python
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
```



## 第一次特异性指标计算(基因组数量)

[tspex: Tissue-specificity calculator (unicamp.br)](https://tspex.lge.ibi.unicamp.br/)



# HMM模型去冗余

## Artificial_protein

```shell
#producing
hmmemit -N 10 -o /home/taoyechen/Balance_process/artificial_protein.fasta /home/taoyechen/Balance_process/hmm_20_combine.hmmre
#aligning
hmmsearch --cpu 180 -E 1e-5 --tblout -o /home/taoyechen/NCBI_Complete_process/test1.txt /home/taoyechen/NCBI_Complete_process/hmm_20_combine.hmm /home/taoyechen/NCBI_Complete_process/artificial_protein.fasta 

nohup hmmsearch --cpu 180 -E 1e-5 --tblout /home/taoyechen/Balance_process/resultrede.txt /home/taoyechen/Balance_process/hmm_20_combine.hmm /home/taoyechen/Balance_process/artificial_protein.fasta >/dev/null 2>&1 &
```

## Result_txt2csv

```python
import pandas as pd

# 读取文本文件内容，跳过注释行
data = []
with open(r'C:\Users\TyCcc\Desktop\test1.txt', 'r') as file:
    for line in file:
        if not line.startswith('#') and line.strip():
            data.append(line.split())

# 定义列名
columns = [
    'seq_id', 'accession', 'hmm_id', 'query_accession',
    'full_E-value', 'full_score', 'full_bias', 'best_1_domain_E-value',
    'best_1_domain_score', 'best_1_domain_bias', 'exp', 'reg', 'clu',
    'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target'
]

# 创建DataFrame
df = pd.DataFrame(data, columns=columns)

# 将DataFrame保存到CSV文件
df.to_csv(r'C:\Users\TyCcc\Desktop\output.csv', index=False)
```

## Cosine_similarity_matrix

```python
import pandas as pd
from itertools import combinations, islice
from scipy.spatial.distance import cosine
import os
from concurrent.futures import ProcessPoolExecutor
import time


# Function to calculate cosine similarity between two vectors
def custom_cosine_similarity(vector1, vector2):
    return 1 - cosine(vector1, vector2)


# Function to process a batch of pairs and calculate cosine similarity
def process_batch(batch, feature_vectors, batch_index):
    results = []
    for pair in batch:
        hmm_id1, hmm_id2 = pair
        # 获取两个hmm_id的seq_id集合
        seq_ids1 = set(feature_vectors[hmm_id1].keys())
        seq_ids2 = set(feature_vectors[hmm_id2].keys())

        # 检查是否存在重复的seq_id
        common_seq_ids = seq_ids1 & seq_ids2
        if common_seq_ids:
            # 如果存在重复的seq_id，计算余弦相似度
            all_seq_ids = seq_ids1 | seq_ids2
            vector1 = [feature_vectors[hmm_id1].get(seq_id, 0) for seq_id in all_seq_ids]
            vector2 = [feature_vectors[hmm_id2].get(seq_id, 0) for seq_id in all_seq_ids]
            similarity = custom_cosine_similarity(vector1, vector2)
        else:
            # 如果不存在重复的seq_id，相似度设为0
            similarity = 0

        results.append((hmm_id1, hmm_id2, similarity))
    # Save results to a temporary CSV file
    temp_file_path = f'temp_results_{batch_index}.csv'
    pd.DataFrame(results, columns=['hmm_id1', 'hmm_id2', 'similarity']).to_csv(temp_file_path, index=False)
    return temp_file_path



def update_matrix_from_file(temp_file_path, cosine_similarity_matrix):
    temp_df = pd.read_csv(temp_file_path)
    # 创建一个空的 DataFrame 用于更新
    update_df = pd.DataFrame(0, index=cosine_similarity_matrix.index, columns=cosine_similarity_matrix.columns, dtype=float)
    for index, row in temp_df.iterrows():
        hmm_id1, hmm_id2, similarity = row['hmm_id1'], row['hmm_id2'], row['similarity']
        if hmm_id1 in cosine_similarity_matrix.index and hmm_id2 in cosine_similarity_matrix.columns:
            update_df.at[hmm_id1, hmm_id2] = similarity
            update_df.at[hmm_id2, hmm_id1] = similarity
    return update_df


# Read the comparison results file into a DataFrame
comparison_results_df = pd.read_csv('/home/taoyechen/Balance_process/out1.csv')

# Filter out rows that do not contain comparison data
comparison_results_df = comparison_results_df[comparison_results_df['target name'].str.startswith('cluster_')]

# Update the score for each sequence in the HMM model's dictionary
feature_vectors = {}
for index, row in comparison_results_df.iterrows():
    hmm_id = row['query name']
    seq_id = row['target name']
    score = row['score']
    if hmm_id not in feature_vectors:
        feature_vectors[hmm_id] = {}
    feature_vectors[hmm_id][seq_id] = score

# Calculate cosine similarity matrix
hmm_ids = list(feature_vectors.keys())
cosine_similarity_matrix = pd.DataFrame(0, index=hmm_ids, columns=hmm_ids, dtype=float)

# Create combinations of HMM ID pairs
pairs_generator = combinations(hmm_ids, 2)

# Divide the pairs into batches
batch_size = 500000  # Adjust the batch size as needed

# Calculate the total number of batches needed
total_pairs = len(hmm_ids) * (len(hmm_ids) - 1) // 2
batches_needed = total_pairs // batch_size + (total_pairs % batch_size > 0)

if __name__ == '__main__':
    start_time = time.perf_counter()

    # Create a generator for batches of data pairs
    batch_generator = (list(islice(pairs_generator, batch_size)) for _ in range(batches_needed))

    try:
        with ProcessPoolExecutor() as executor:
            # Process batches and save results to temporary files
            temp_file_paths = list(executor.map(process_batch, batch_generator, [feature_vectors] * batches_needed, range(batches_needed)))

            # Submit all file processing tasks
            futures = [executor.submit(update_matrix_from_file, temp_file_path, cosine_similarity_matrix) for temp_file_path in temp_file_paths]

            # Initialize an empty DataFrame for cumulative updates
            cumulative_update_df = pd.DataFrame(0, index=cosine_similarity_matrix.index, columns=cosine_similarity_matrix.columns, dtype=float)
            for future in futures:
                # Wait for each task to complete and update the cumulative DataFrame
                cumulative_update_df += future.result()

    finally:
        # Delete temporary files
        for temp_file_path in temp_file_paths:
            if os.path.exists(temp_file_path):
                os.remove(temp_file_path)

    # Update the main DataFrame
    cosine_similarity_matrix.update(cumulative_update_df)

    # Save the final DataFrame to a CSV file
    cosine_similarity_matrix.to_csv('/home/taoyechen/Balance_process/cosine_similarity_matrix.csv')
    print("Cosine similarity matrix has been calculated and saved.")
    end_time = time.perf_counter()
    run_time_minutes = (end_time - start_time) / 60
    print(f"运行时间：{run_time_minutes:.2f}分钟")
```

## Sparse_matrix

```python
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import igraph as ig
import leidenalg

# 读取CSV文件
df = pd.read_csv('/home/taoyechen/parameter/allcosine_similarity_matrix.csv', index_col=0)

# 获取hmmid列表
hmmids = df.index.tolist()

# 将对角线上的值设置为1
np.fill_diagonal(df.values, 1)

# 创建CSR格式的稀疏矩阵
sparse_matrix = csr_matrix(df.values)

# 使用igraph创建图
sources, targets = sparse_matrix.nonzero()
weights = sparse_matrix[sources, targets].tolist()[0]
g = ig.Graph(directed=False)
g.add_vertices(hmmids)
g.add_edges(list(zip(sources.tolist(), targets.tolist())))
g.es['weight'] = weights

# 设置分辨率参数
resolution_parameter = 0.25  # 可以根据需要调整这个值

# 应用Leiden算法进行社区检测，包含分辨率参数
partition = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, resolution_parameter=resolution_parameter, weights='weight')

# 创建一个DataFrame来存储hmmid和它们的社区编号
community_df = pd.DataFrame({'hmmid': hmmids, 'community': partition.membership})

# 将DataFrame保存为CSV文件
community_df.to_csv('/home/taoyechen/parameter/leiden_clustering_results.csv', index=False)

# 打印操作成功的消息
print("聚类结果已保存到 'leiden_clustering_results.csv' 文件中。")

```

Leiden聚类算法中的分辨率参数是用来控制社区（即聚类）检测的粒度的。这个参数直接影响到算法划分社区的大小和数量。具体来说：

- **较高的分辨率**：会导致算法识别出更多的社区，但每个社区的规模较小。这意味着算法会将数据集划分成更细致的组，每个组包含的节点较少。
- **较低的分辨率**：则会导致算法识别出较少的社区，但每个社区的规模较大。在这种情况下，算法倾向于将节点聚集成较大的组。

  分辨率参数的作用可以类比于一个阈值，它决定了社区内部的连接密度应该至少达到多少，而社区之间的连接密度应该低于多少。因此，通过调整分辨率参数，可以控制聚类的粗细程度，从而影响最终的聚类结果。

在实际应用中，选择合适的分辨率参数需要根据数据的特性和分析的目的来决定。

## Hmm_cluster_select

聚类后在每簇中选择香浓特异性指标值(使用第一次计算得到的特异性指标)最大的hmm模型,得到hmmid.txt

```python
import pandas as pd

# 读取CSV文件
df = pd.read_csv(r"C:\Users\TyCcc\Desktop\leiden_clustering_results.csv")

# 按community分组，并找到每个分组中shannon_specificity的最大值
max_shannon = df.groupby('community')['shannon_specificity'].transform(max)

# 创建一个布尔掩码，标记每个分组中shannon_specificity等于最大值的行
mask = df['shannon_specificity'] == max_shannon

# 在每个community分组中，如果有多个最大值，则随机选择一行
# 我们可以通过对满足条件的行进行分组并随机选择一行来实现
# 使用cumcount为每个分组中的行分配一个唯一的序号，并基于序号进行筛选
# 这里，我们假设如果有多行具有相同的最大值，我们只需要其中一行
grouped = df[mask].groupby('community', as_index=False, group_keys=False)
random_rows = []
for name, group in grouped:
    if len(group) > 1:
        # 如果有多行具有相同的最大值，随机选择一行
        random_row = group.sample(n=1)
        random_rows.append(random_row)
    else:
        # 如果只有一行具有最大值，直接添加
        random_rows.append(group)

    # 将列表中的DataFrame合并成一个DataFrame
result_df = pd.concat(random_rows, ignore_index=True)

# 将结果保存到新的CSV文件中
result_df.to_csv(r"C:\Users\TyCcc\Desktop\balance_hmm_dere.csv", index=False)

print("筛选后的数据已保存到filtered_data.csv文件中。")
```



## Hmm_clustered_link

```shell
#!/bin/bash  
  
# 设置目录和文件路径  
HMMID_FILE="/home/taoyechen/Balance_process/hmmid.txt"  
HMM_MODELS_DIR="/home/taoyechen/Balance_process/hmm_20"  
LINKS_DIR="/home/taoyechen/Balance_process/hmm_clustered"  
  
# 确保链接目录存在  
mkdir -p "$LINKS_DIR"  
  
# 读取hmmid.txt中的每一行  
while IFS= read -r hmmid; do  
    # 构建hmm模型文件的完整路径  
    HMM_FILE="${HMM_MODELS_DIR}${hmmid}.hmm"  
    # 构建链接的完整路径  
    LINK_FILE="${LINKS_DIR}${hmmid}.hmm"  
      
    # 检查hmm模型文件是否存在  
    if [ -f "$HMM_FILE" ]; then  
        # 如果链接文件已存在，则删除它（如果需要的话）  
        if [ -e "$LINK_FILE" ]; then  
            rm "$LINK_FILE"  
        fi  
        # 创建到hmm模型文件的符号链接  
        ln -s "$HMM_FILE" "$LINK_FILE"  
        echo "Created link: $LINK_FILE -> $HMM_FILE"  
    else  
        echo "No hmm file found for $hmmid"  
    fi  
done < "$HMMID_FILE"  
  
echo "All links have been created or checked."
```



# Marker获取

## Protein_num(domz)

```python
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

```



## 第二次特异性指标计算(蛋白数量)

[tspex: Tissue-specificity calculator (unicamp.br)](https://tspex.lge.ibi.unicamp.br/)



## Marker_select

```python
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

```



## Hmm_link_Marker

```python
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


```



# Feature_get

## Marker_hmm_align

```shell
#!/bin/bash

PATH_HMM="/home/taoyechen/classification_pre/New_feature_get/pm.hmm"
NONPATH_HMM="/home/taoyechen/classification_pre/New_feature_get/nm.hmm"
PATH_SEQ="/home/taoyechen/classification_pre/hmm_align/path/prodigal_protein"
NONPATH_SEQ="/home/taoyechen/classification_pre/hmm_align/non_p/prodigal_protein"

P2P="/home/taoyechen/classification_pre/New_feature_get/p_2_pm"
P2N="/home/taoyechen/classification_pre/New_feature_get/p_2_nm"
N2P="/home/taoyechen/classification_pre/New_feature_get/n_2_pm"
N2N="/home/taoyechen/classification_pre/New_feature_get/n_2_nm"

# 定义比对函数
run_hmmsearch() {
  local hmm_model=$1
  local sequence_dir=$2
  local output_dir=$3

  # 确保输出目录存在
  mkdir -p "$output_dir"

  # 遍历文件夹中的所有序列文件
  for sequence_file in "$sequence_dir"/*
  do
    # 获取序列文件的基本名
    local base_name=$(basename "$sequence_file")
    
    # 在后台运行hmmsearch并将结果输出到指定路径
    hmmsearch --cpu 2 -E 1e-5 -o "${output_dir}/${base_name}.txt" "$hmm_model" "$sequence_file" &
  done
  # 等待所有后台进程完成
  wait
}

# 并行运行比对函数
run_hmmsearch "$PATH_HMM" "$PATH_SEQ" "$P2P" &
run_hmmsearch "$PATH_HMM" "$NONPATH_SEQ" "$N2P" &
run_hmmsearch "$NONPATH_HMM" "$NONPATH_SEQ" "$N2N" &
run_hmmsearch "$NONPATH_HMM" "$PATH_SEQ" "$P2N" &
wait

echo "所有比对已完成。"
```

## align_result

```python
import re
import csv
import os
import glob
import pandas as pd


def extract_info(filenames, csv_file_path):
    # 初始化字典来存储每个query_id的domZ值列表
    query_info = {}

    # 循环处理每个文件
    for filename in filenames:
        with open(filename, 'r') as file:
            data = file.read()

        # 提取所有query_id和domZ值
        query_ids = re.findall(r'Query:\s+(\S+)', data)
        domZ_values = re.findall(r'Domain search space\s+\(domZ\):\s+(\d+)', data)

        # 将每个query_id的domZ值添加到字典中
        for query_id, domZ_value in zip(query_ids, domZ_values):
            if query_id not in query_info:
                query_info[query_id] = []
            query_info[query_id].append(domZ_value)

    # 提取文件名的指定部分作为列标题
    column_titles = []
    for f in filenames:
        match = re.search(r'(.*?_.*?)_', os.path.basename(f))
        if match:
            column_titles.append(match.group(1))
        else:
            column_titles.append(os.path.basename(f))  # 如果没有匹配，使用整个文件名

    # 打开CSV文件并准备写入
    with open(csv_file_path, 'w', newline='') as csvfile:
        # 字段名为'Query ID'和每个文件的指定部分
        fieldnames = ['Query ID'] + column_titles
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        # 为每个query_id写入一行，包含所有domZ值
        for query_id, domZ_values in query_info.items():
            row = {'Query ID': query_id}
            # 确保每个domZ值都在正确的列中
            for i, domZ_value in enumerate(domZ_values):
                row[fieldnames[i+1]] = domZ_value
            writer.writerow(row)
# 文件夹路径
folder_path1 = r'/home/taoyechen/classification_pre/New_feature_get/p_2_pm'

folder_path2 = r'/home/taoyechen/classification_pre/New_feature_get/n_2_pm'

folder_path3 = r'/home/taoyechen/classification_pre/New_feature_get/n_2_nm'

folder_path4 = r'/home/taoyechen/classification_pre/New_feature_get/p_2_nm'

csv_file_path1 = r'/home/taoyechen/classification_pre/New_feature_get/p_2_pm.csv'
csv_file_path2 = r'/home/taoyechen/classification_pre/New_feature_get/n_2_pm.csv'
csv_file_path3 = r'/home/taoyechen/classification_pre/New_feature_get/n_2_nm.csv'
csv_file_path4 = r'/home/taoyechen/classification_pre/New_feature_get/p_2_nm.csv'

# 获取文件夹中所有.txt文件的列表
filenames1 = glob.glob(os.path.join(folder_path1, '*.txt'))
filenames2 = glob.glob(os.path.join(folder_path2, '*.txt'))
filenames3 = glob.glob(os.path.join(folder_path3, '*.txt'))
filenames4 = glob.glob(os.path.join(folder_path4, '*.txt'))
# 使用函数
extract_info(filenames1, csv_file_path1)
extract_info(filenames2, csv_file_path2)
extract_info(filenames3, csv_file_path3)
extract_info(filenames4, csv_file_path4)




def merge_csv_columns(csv_file1, csv_file2, insert_after_column_index, output_file_name, output_directory):
    # 读取两个CSV文件
    df1 = pd.read_csv(csv_file1)
    df2 = pd.read_csv(csv_file2)

    # 获取csv_file2中要插入的列
    columns_to_insert = df2.iloc[:, insert_after_column_index:insert_after_column_index+2]

    # 在df1中创建新的列名，以便在插入时不会覆盖原有的列
    new_column_names = [f'{col}' for col in columns_to_insert.columns]

    # 重命名df2的列，以便在合并时保持它们的唯一性
    columns_to_insert.columns = new_column_names

    # 将新列插入到df1的指定位置
    for i, col in enumerate(new_column_names, start=insert_after_column_index):
        df1.insert(i, col, columns_to_insert[col])

    # 构建输出文件的完整路径
    output_file_path = f'{output_directory}\\{output_file_name}'

    # 将合并后的DataFrame保存到指定的输出文件
    df1.to_csv(output_file_path, index=False)

    return output_file_path

# 调用函数
output_path = merge_csv_columns(
    csv_file1=csv_file_path1,
    csv_file2=r'D:\workplace_pycharm\pythonProject1\forest\pathogen_specific_or.csv',
    insert_after_column_index=1,
    output_file_name='pps.csv',
    output_directory=r'C:\Users\TyCcc\Desktop'
)

output_path = merge_csv_columns(
    csv_file1=csv_file_path2,
    csv_file2=r'D:\workplace_pycharm\pythonProject1\forest\pathogen_specific_or.csv',
    insert_after_column_index=1,
    output_file_name='nps.csv',
    output_directory=r'C:\Users\TyCcc\Desktop'
)

output_path = merge_csv_columns(
    csv_file1=csv_file_path3,
    csv_file2=r'D:\workplace_pycharm\pythonProject1\forest\non_pathogen_specific_or.csv',
    insert_after_column_index=1,
    output_file_name='nns.csv',
    output_directory=r'C:\Users\TyCcc\Desktop'
)

output_path = merge_csv_columns(
    csv_file1=csv_file_path4,
    csv_file2=r'D:\workplace_pycharm\pythonProject1\forest\non_pathogen_specific_or.csv',
    insert_after_column_index=1,
    output_file_name='pns.csv',
    output_directory=r'C:\Users\TyCcc\Desktop'
)

print(f'合并后的文件已保存到: {output_path}')
```



```python
import os  
import re  
import csv  
  
# 假设这是包含hmmid的txt文件路径  
hmmid_file_path = '/home/taoyechen/Balance_process/hmmid.txt'  
# 假设这是包含比对结果文件夹的父目录路径  
alignment_folder_path = '/home/taoyechen/Balance_process/hp_result'  
# CSV文件输出路径  
csv_output_path = '/home/taoyechen/Balance_process/hp_result_domZ.csv'  
  
# 读取hmmid文件并获取所有的hmmid  
with open(hmmid_file_path, 'r') as hmmid_file:  
    hmmids = [line.strip() for line in hmmid_file.readlines()]  
  
# 初始化一个空列表来存储所有需要处理的文件路径  
filenames = []  
  
# 遍历所有的hmmid  
for hmmid in hmmids:  
    # 构造hmmid对应的文件夹路径  
    hmmid_folder_path = os.path.join(alignment_folder_path, hmmid)  
      
    # 检查文件夹是否存在  
    if os.path.exists(hmmid_folder_path) and os.path.isdir(hmmid_folder_path):  
        # 遍历文件夹中的txt文件（比对结果文件）  
        for alignment_file in os.listdir(hmmid_folder_path):  
            if alignment_file.endswith('.txt'):  
                filenames.append(os.path.join(hmmid_folder_path, alignment_file))  
  
# 调用函数处理文件并写入CSV  
extract_info(filenames, csv_output_path)  
  
# extract_info函数的定义（保持与之前相同）  
def extract_info(filenames, csv_file_path):  
    # 初始化字典来存储每个query_id的domZ值列表  
    query_info = {}  
  
    # 循环处理每个文件  
    for filename in filenames:  
        with open(filename, 'r') as file:  
            data = file.read()  
  
        # 提取所有query_id和domZ值  
        query_ids = re.findall(r'Query:\s+(\S+)', data)  
        domZ_values = re.findall(r'Domain search space\s+\(domZ\):\s+(\d+)', data)  
  
        # 将每个query_id的domZ值添加到字典中  
        for query_id, domZ_value in zip(query_ids, domZ_values):  
            if query_id not in query_info:  
                query_info[query_id] = []  
            query_info[query_id].append(domZ_value)  
  
    # 提取文件名的指定部分作为列标题  
    column_titles = []  
    for f in filenames:  
        match = re.search(r'(.*?_.*?)_', os.path.basename(f))  
        if match:  
            column_titles.append(match.group(1))  
        else:  
            column_titles.append(os.path.basename(f))  # 如果没有匹配，使用整个文件名  
  
    # 打开CSV文件并准备写入  
    with open(csv_file_path, 'w', newline='') as csvfile:  
        # 字段名为'Query ID'和每个文件的指定部分  
        fieldnames = ['Query ID'] + column_titles  
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)  
  
        writer.writeheader()  
        # 为每个query_id写入一行，包含所有domZ值  
        for query_id, domZ_values in query_info.items():  
            row = {'Query ID': query_id}  
            # 确保每个domZ值都在正确的列中  
            for i, domZ_value in enumerate(domZ_values):  
                row[fieldnames[i+1]] = domZ_value  
            writer.writerow(row)  
  
# 脚本执行完毕
```

## Features_compute

### All

```python
import pandas as pd
import numpy as np

# 加载CSV文件
df = pd.read_csv(r"C:\Users\TyCcc\Desktop\hp_result_M.csv")

# 初始化一个空的DataFrame来保存结果
result_df = pd.DataFrame(index=df.columns[3:].tolist(),
                         columns=['NHP_Count', 'HP_Count', 'NHP_SPM_Median', 'HP_SPM_Median', 'NHP_SPM_Total',
                                  'HP_SPM_Total', 'Sigmoid_NHP_over_HP', 'Sigmoid_HP_over_NHP'])


# 定义sigmoid函数
def sigmoid(x):
    return 1 / (1 + np.exp(-x))


# 遍历每个基因组列
for genome_column in df.columns[3:]:
    nhp_counts = 0
    hp_counts = 0
    nhp_spm_values = []
    hp_spm_values = []
    nhp_spm_total = 0
    hp_spm_total = 0

    # 遍历每一行
    for index, row in df.iterrows():
        if row[genome_column] > 0:  # 检测到hmm模型
            if row['marker'] == 0:  # NHP标签
                nhp_counts += row[genome_column]
                nhp_spm_values.append(row['NHP_SPM'])
                nhp_spm_total += row['NHP_SPM']
            elif row['marker'] == 1:  # HP标签
                hp_counts += row[genome_column]
                hp_spm_values.append(row['HP_SPM'])
                hp_spm_total += row['HP_SPM']

                # 计算中位数（如果列表不为空）
    nhp_spm_median = np.median(nhp_spm_values) if nhp_spm_values else np.nan
    hp_spm_median = np.median(hp_spm_values) if hp_spm_values else np.nan

    # 计算sigmoid值
    sigmoid_nhp_over_hp = sigmoid(nhp_spm_total - hp_spm_total) if not np.isnan(nhp_spm_total) and not np.isnan(
        hp_spm_total) else np.nan
    sigmoid_hp_over_nhp = sigmoid(hp_spm_total - nhp_spm_total) if not np.isnan(hp_spm_total) and not np.isnan(
        nhp_spm_total) else np.nan

    # 将结果保存到DataFrame中
    result_df.at[genome_column, 'NHP_Count'] = nhp_counts
    result_df.at[genome_column, 'HP_Count'] = hp_counts
    result_df.at[genome_column, 'NHP_SPM_Median'] = nhp_spm_median
    result_df.at[genome_column, 'HP_SPM_Median'] = hp_spm_median
    result_df.at[genome_column, 'NHP_SPM_Total'] = nhp_spm_total
    result_df.at[genome_column, 'HP_SPM_Total'] = hp_spm_total
    result_df.at[genome_column, 'Sigmoid_NHP_over_HP'] = sigmoid_nhp_over_hp
    result_df.at[genome_column, 'Sigmoid_HP_over_NHP'] = sigmoid_hp_over_nhp

# 将结果保存到新的CSV文件中
result_df.to_csv(r"C:\Users\TyCcc\Desktop\hp_feature.csv")

```

### Features_pre

```python
import pandas as pd

# 读取CSV文件
df = pd.read_csv(r'C:\Users\TyCcc\Desktop\n_2_pm.csv')

# 准备一个空的DataFrame来存储结果
results = pd.DataFrame()

# 计算每个序列的中位数和总和
for sequence in df.columns[2:]:
    # 筛选出数量大于0的marker
    filtered_df = df[df[sequence] > 0]
    # 如果filtered_df不为空，计算特异性指标值的中位数；否则，中位数为0
    median_value = filtered_df['pathogen_SPM'].median() if not filtered_df.empty else 0
    # 计算marker数量的总和
    marker_count = df[sequence].sum()
    # 计算特异性指标值的总和
    specificity_sum = filtered_df['pathogen_SPM'].sum() if not filtered_df.empty else 0
    # 将结果添加到results DataFrame中
    results = results._append({'sequence': sequence, 'p_median_value': median_value, 'p_marker_count': marker_count, 'p_specificity_sum': specificity_sum}, ignore_index=True)

# 将结果输出到CSV文件
results.to_csv(r'C:\Users\TyCcc\Desktop\n2p.csv', index=False)
```

### Feature_compute

```python
import pandas as pd
import numpy as np

# 定义sigmoid函数
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# 读取CSV文件
df = pd.read_csv(r'C:\Users\TyCcc\Desktop\features.csv')

# 计算p_marker_freq和n_marker_freq
df['p_marker_freq'] = df['p_marker_count'] / df['seq_num']
df['n_marker_freq'] = df['n_marker_count'] / df['seq_num']

# 计算p_vs_n_score_logistic和n_vs_p_score_logistic
df['p_vs_n_score_logistic'] = sigmoid(df['p_specificity_sum'] - df['n_specificity_sum'])
df['n_vs_p_score_logistic'] = sigmoid(df['n_specificity_sum'] - df['p_specificity_sum'])

# 将结果保存到新的CSV文件
df.to_csv('updated_file.csv', index=False)

```

## seq_num(pre)

```python
import re
import csv
import os
import glob

def extract_info(filenames):
    # 初始化字典来存储每个query_id的Target sequences值列表
    query_info = {}

    # 循环处理每个文件
    for filename in filenames:
        with open(filename, 'r') as file:
            data = file.read()

        # 提取所有query_id和Target sequences值
        query_ids = re.findall(r'Query:\s+(\S+)', data)
        target_sequences_values = re.findall(r'Target sequences:\s+(\d+)', data)

        # 将每个query_id的Target sequences值添加到字典中
        for query_id, target_sequences_value in zip(query_ids, target_sequences_values):
            if query_id not in query_info:
                query_info[query_id] = []
            query_info[query_id].append(target_sequences_value)

    # 提取文件名的指定部分作为列标题
    column_titles = []
    for f in filenames:
        match = re.search(r'(.*?_.*?)_', os.path.basename(f))
        if match:
            column_titles.append(match.group(1))
        else:
            column_titles.append(os.path.basename(f))  # 如果没有匹配，使用整个文件名

    # 打开CSV文件并准备写入
    with open('p_num.csv', 'w', newline='') as csvfile:
        # 字段名为'Query ID'和每个文件的指定部分
        fieldnames = ['Query ID'] + column_titles
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        # 为每个query_id写入一行，包含所有Target sequences值
        for query_id, target_sequences_values in query_info.items():
            row = {'Query ID': query_id}
            # 确保每个Target sequences值都在正确的列中
            for i, target_sequences_value in enumerate(target_sequences_values):
                row[fieldnames[i+1]] = target_sequences_value
            writer.writerow(row)

# 文件夹路径
folder_path = r'/home/taoyechen/classification_pre/feature_get/result/p_2_nm'
# 获取文件夹中所有.txt文件的列表
filenames = glob.glob(os.path.join(folder_path, '*.txt'))
# 使用函数
extract_info(filenames)

```

# 模型训练

## XGBOOST

```python
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.metrics import f1_score, precision_recall_curve, auc, roc_curve, confusion_matrix, matthews_corrcoef
import matplotlib.pyplot as plt
from sklearn.metrics import recall_score  # Import missing function


# 读取数据集
dataset_path = r'C:\Users\TyCcc\Desktop\feature_per.csv' # 替换为你的数据集路径
df = pd.read_csv(dataset_path)

# 分割特征和目标变量
X = df.drop(columns=['label', 'sequence'])  # 假设'seq_name'是不需要的列
y = df['label']

# 数据标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# 定义要搜索的参数网格
param_grid = {
    'n_estimators': [100, 200],
    'max_depth': [None, 10, 20, 30, 40],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
}

# 创建随机森林分类器实例
rf_classifier = RandomForestClassifier(random_state=42)

# 使用GridSearchCV进行参数搜索
grid_search = GridSearchCV( n_jobs=-1, estimator=rf_classifier, param_grid=param_grid, cv=5, scoring='accuracy')
grid_search.fit(X_train, y_train)

# 获取最佳参数
best_params = grid_search.best_params_
print("Best Parameters: ", best_params)

# 使用最佳参数重新训练模型
best_clf = grid_search.best_estimator_

# 在测试集上进行预测
y_pred = best_clf.predict(X_test)

# 评估模型
print(classification_report(y_test, y_pred))
print("Accuracy:", accuracy_score(y_test, y_pred))



# 计算F1-macro
f1_macro = f1_score(y_test, y_pred, average='macro')
print("F1-macro:", f1_macro)

# 计算PR-AUC
precision, recall, thresholds = precision_recall_curve(y_test, best_clf.predict_proba(X_test)[:, 1])
pr_auc = auc(recall, precision)
print("PR-AUC:", pr_auc)

# 计算ROC-AUC
fpr, tpr, roc_thresholds = roc_curve(y_test, best_clf.predict_proba(X_test)[:, 1])
roc_auc = auc(fpr, tpr)
print("ROC-AUC:", roc_auc)

# 计算Sensitivity (Recall)
sensitivity = recall_score(y_test, y_pred)
print("Sensitivity:", sensitivity)

# 计算Specificity
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
specificity = tn / (tn + fp)
print("Specificity:", specificity)

# 计算MCC
mcc = matthews_corrcoef(y_test, y_pred)
print("MCC:", mcc)

# 绘制ROC曲线
plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()

# 绘制PR曲线
plt.figure()
plt.plot(recall, precision, color='blue', lw=2, label='Precision-Recall curve (area = %0.2f)' % pr_auc)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend(loc="lower left")
plt.show()
```

## compare

```python
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
import xgboost as xgb
import numpy as np
from sklearn.ensemble import VotingClassifier


# 加载数据
data = pd.read_csv(r'/home/taoyechen/Balance_process/all_feature_592.csv')

# 分离特征和标签
X = data.drop(['Genome', 'label'], axis=1)  # 特征列
y = data['label']  # 标签列


# 数据分割：分层采样
X_train_stratified, X_temp, y_train_stratified, y_temp = train_test_split(X, y, test_size=0.4, random_state=42,
                                                                          stratify=y)
X_val_stratified, X_test_stratified, y_val_stratified, y_test_stratified = train_test_split(X_temp, y_temp,
                                                                                            test_size=0.5,
                                                                                            random_state=43,
                                                                                            stratify=y_temp)

# 数据分割：随机采样
X_train_random, X_temp, y_train_random, y_temp = train_test_split(X, y, test_size=0.4, random_state=44)
X_val_random, X_test_random, y_val_random, y_test_random = train_test_split(X_temp, y_temp, test_size=0.5,
                                                                            random_state=45)

# 定义模型
models = {
    'LogisticRegression': LogisticRegression(),
    'RandomForest': RandomForestClassifier(),
    'SVM': SVC(),
    'XGBoost': xgb.XGBClassifier(use_label_encoder=False, eval_metric='mlogloss')
}

# 定义参数网格
param_grids = {
    'LogisticRegression': {
        'C': [0.01, 0.1, 1, 10, 100],
        'penalty': ['l1', 'l2'],
        'solver': ['liblinear', 'saga']
    },
    'RandomForest': {
        'n_estimators': [100, 200, 300],
        'max_depth': [None, 10, 20, 30],
        'min_samples_split': [2, 5, 10]
    },
    'SVM': {
        'C': [0.1, 1, 10, 100],
        'gamma': ['scale', 'auto', 0.01, 0.001, 0.0001],
        'kernel': ['rbf', 'linear', 'poly', 'sigmoid']
    },
    'XGBoost': {
        'n_estimators': [100, 200, 300],
        'learning_rate': [0.01, 0.05, 0.1],
        'max_depth': [3, 5, 7],
        'colsample_bytree': [0.5, 0.7],
        'gamma': [0, 0.25, 1.0]
    }
}


# 训练并调优模型的函数
def train_and_evaluate_model(model_name, model, param_grid, X_train, y_train, X_test, y_test):
    grid_search = GridSearchCV(model, param_grid, cv=5, scoring='accuracy', verbose=1)
    grid_search.fit(X_train, y_train)
    best_model = grid_search.best_estimator_
    y_pred = best_model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Best parameters for {model_name}: {grid_search.best_params_}")
    print(f"Test Accuracy for {model_name}: {accuracy}")
    return best_model, accuracy


# 初始化字典来存储最佳模型和它们的准确率
best_models = {}
accuracies = {}

# 使用分层采样数据训练并评估模型
for model_name, model in models.items():
    param_grid = param_grids[model_name]
    best_model, accuracy = train_and_evaluate_model(model_name, model, param_grid,
                                                    X_train_stratified, y_train_stratified,
                                                    X_test_stratified, y_test_stratified)
    best_models[f'{model_name}_stratified'] = best_model
    accuracies[f'{model_name}_stratified_accuracy'] = accuracy

# 使用随机采样数据训练并评估模型
for model_name, model in models.items():
    param_grid = param_grids[model_name]
    best_model, accuracy = train_and_evaluate_model(model_name, model, param_grid,
                                                    X_train_random, y_train_random,
                                                    X_test_random, y_test_random)
    best_models[f'{model_name}_random'] = best_model
    accuracies[f'{model_name}_random_accuracy'] = accuracy

# 打印所有模型的准确率
print("\nModel Accuracies with Stratified Sampling:")
for model, accuracy in accuracies.items():
    if 'stratified' in model:
        print(f"{model}: {accuracy}")

print("\nModel Accuracies with Random Sampling:")
for model, accuracy in accuracies.items():
    if 'random' in model:
        print(f"{model}: {accuracy}")



# 假设您已经有了分别训练的四个模型，存储在best_models字典中
# best_models = {
#     'LogisticRegression_stratified': logistic_regression_model,
#     'RandomForest_stratified': random_forest_model,
#     'SVM_stratified': svm_model,
#     'XGBoost_stratified': xgboost_model
# }

# 提取要集成的模型
logistic_regression_model = best_models['LogisticRegression_stratified']
random_forest_model = best_models['RandomForest_stratified']
svm_model = best_models['SVM_stratified']
xgboost_model = best_models['XGBoost_stratified']

# 创建集成分类器
voting_classifier = VotingClassifier(
    estimators=[
        ('log_reg', logistic_regression_model),
        ('rand_forest', random_forest_model),
        ('svm', svm_model),
        ('xgboost', xgboost_model)
    ],
    voting='soft'  # 或者使用 'hard' 投票，根据您的需求选择
)

# 使用之前分割出来的验证集来评估集成模型的性能
voting_classifier.fit(X_train_stratified, y_train_stratified)  # 假设您使用的是分层采样数据
y_pred = voting_classifier.predict(X_val_stratified)
accuracy = accuracy_score(y_val_stratified, y_pred)

print(f"Validation Accuracy for Voting Classifier: {accuracy}")

# 如果需要，您还可以在测试集上评估模型性能
y_test_pred = voting_classifier.predict(X_test_stratified)
test_accuracy = accuracy_score(y_test_stratified, y_test_pred)
print(f"Test Accuracy for Voting Classifier: {test_accuracy}")

```



```shell
cut -f2 /home/shiyuyu/24_0425/calssification_pre/dere_new/drclu.tsv | sed 's/^/>/' > /dev/shm/ID && sed '/^>/ {s/>lcl|/>/;s/ .*//}' /home/shiyuyu/24_0425/calssification_pre/combined_proteins.fasta | awk 'BEGIN{RS=">";FS="\n";ORS=""} NR==FNR{a[$1];next} $1 in a{print ">"$0}' /dev/shm/ID - > /dev/shm/extracted.fasta
```

1. 获取数据：
   - 数据1：从NCBI获取的ASSEMBLY LEVEL为complete的细菌基因组元数据。
   - 数据2：从GCM获取的ASSEMBLY LEVEL为complete的人类病原菌基因组元数据。
   - 数据3：从BV-BRC和IMG获取的ASSEMBLY LEVEL为complete的宿主为人类的细菌基因组元数据。
   - 数据4：中国疾控中心发布的《人间传染的病原微生物目录》中病原菌的属种信息。
2. 筛选主要病原菌列表信息：
   1. 将数据2中的元数据作为主要病原菌列表。
   2. 使用Genbank Accession将数据1与数据2进行比较，取其重叠部分，作为主要病原菌列表信息。
3. 筛选补充病原菌列表信息：
   1. 将数据3中两个数据库元数据根据已有文章的元数据注释筛选方法进行筛选。
   2. 将筛选出的数据与数据4中的种属信息进行比较，取重叠部分，作为补充病原菌列表信息。
4. 整合并去重：
   1. 将步骤2中获取的主要病原菌列表信息与步骤3中获取的补充病原菌列表信息进行整合并去重，形成数据5。
   2. 利用Genbank Accession将数据1与数据5进行比较，取其重叠部分，作为最终病原菌列表信息（数据6）。
5. 筛选非病原菌列表信息：
   1. 使用Genbank Accession将数据1与数据6进行比较，取数据1中的反集，作为非病原菌列表信息。
   2. 将该筛选出的反集数据与数据4中的种属信息进行比较，取不重叠部分，作为最终的非病原菌列表信息。

以下是基于上述步骤的流程图：
                                  +----------------------+
                                 | 数据1: NCBI元数据     |
                                 +----------+-----------+
                                            |
                                            v
                                 +----------+-----------+
                                 | 数据2: GCM元数据     |
                                 +----------+-----------+
                                            |
                                            v
                                 +----------+-----------+
                                 | 主要病原菌列表信息     |
                                 +----------+-----------+
                                            |
                                            v
                                 +----------+-----------+
                                 | 数据3: BV-BRC和IMG元数据|
                                 +----------+-----------+
                                            |
                                            v
                                 +----------+-----------+
                                 | 补充病原菌列表信息     |
                                 +----------+-----------+
                                            |
                                            v
                                 +----------+-----------+
                                 | 整合并去重后的病原菌列表|
                                 +----------+-----------+
                                            |
                                            v
                                 +----------+-----------+
                                 | 最终病原菌列表 (数据6)|
                                 +----------+-----------+
                                            |
                                            v
                                 +----------+-----------+
                                 | 反集 (数据1 - 数据6)  |
                                 +----------+-----------+
                                            |
                                            v
                                 +----------+-----------+
                                 | 最终非病原菌列表      |
                                 +----------------------+

使用的数据为： 

1.National Center for Biotechnology Information获取的ASSEMBLY LEVEL为complete的细菌基因组元数据 

2.Global Catalogue of Microorganisms获取的ASSEMBLY LEVEL为complete的人类病原菌基因组元数据 

3.BACTERIAL AND VIRAL BIOINFORMATICS RESOURCE CENTER、Integrated Microbial Genomes&Microbiomes获取的ASSEMBLY LEVEL为complete的宿主为人类的细菌基因组元数据 

4.中国疾控中心发布的《人间传染的病原微生物目录》中病原菌的属种信息 

筛选步骤为： 

1.将数据2中的元数据作为主要病原菌列表信息，利用Genbank Accession将数据1与数据2元信息进行比较，取其重叠部分作为主要病原菌列表信息。 

2.将数据3中两个数据库元数据按照已有文章的元数据注释筛选方法进行筛选，并将该筛选出的数据与数据4中的种属信息进行比较，取重叠部分作为补充病原菌列表信息。 

3.将步骤1与步骤2获取到的主要病原菌列表信息与补充病原菌列表信息进行整合并去重作为数据5，利用Genbank Accession将数据1与数据5元信息进行比较，取其重叠部分作为最终病原菌列表信息数据6。

 4.利用Genbank Accession将数据1与数据6元信息进行比较，取数据1中的反集作为非病原菌列表信息，并将该筛选出的数据与数据4中的种属信息进行比较，取不重叠部分作为非病原菌列表信息。 

经过该四个步骤，最终得到病原菌列表和非病原菌列表。
 简化的步骤图表示你的过程：

1. **数据准备**

o  病原菌列表

o  非病原菌列表

o  筛选基因组大小为中位数的蛋白质数据

1. **去冗余与聚类**

o  使用mmseqs2的linclust功能去冗余

o  使用mmseqs2的cluster功能聚类

1. **模型构建**

o  筛选序列数目超过20条的簇

o  使用kalign对齐蛋白质序列

o  使用hmmbuild构建hmm模型

1. **特异性指标计算**

o  hmm模型与病原菌/非病原菌蛋白比对

o  使用tspex计算特异性指标

1. **模型去冗余**

o  使用hmmemit模拟人工蛋白质

o  人工蛋白质与hmm模型比对

o  计算余弦相似度

o  使用Leiden算法聚类hmm模型

o  根据特异性指标筛选去冗余模型

1. **特征选择与模型评估**

o  使用去冗余后的hmm模型与蛋白数据比对

o  计算特异性指标作为marker筛选阈值

o  筛选病原菌/非病原菌的蛋白marker

o  计算特征值用于机器学习

o  使用机器学习模型进行分类器训练

使用的数据为：从病原菌列表和非病原菌列表中依据属种信息筛选出基因组大小为该种中中位数的蛋白质数据

步骤为：

1.将蛋白质序列使用mmseqs2软件的linclust功能进行去冗余操作

2.将去冗余后的蛋白质序列使用mmseqs2软件的cluster功能进行聚类操作

3.使用聚类后的蛋白质序列簇筛选出序列数目超过20条的簇，将筛选后的簇使用kalign软件进行对齐，对齐后使用hmm软件的hmmbulid功能将这些蛋白质簇转化为hmm模型

4.将得到的hmm模型与病原菌和非病原菌的蛋白数据进行比对，并根据每个模型是否比对上该蛋白数据，使用tspex软件计算每个hmm模型对于病原菌和非病原菌的特异性指标

5.为了避免出现hmm模型中出现多个簇实际为一簇的情况，对hmm模型进行去冗余操作：

 （1）.使用hmm软件hmmemit功能依据hmm模型模拟出人工蛋白质

 （2）.将模拟出的人工蛋白质与hmm模型进行比对，根据比对结果计算出每个hmm模型的余弦相似度

 （3）.利用余弦相似度矩阵使用Leiden算法进行聚类操作

 （4）.取每个hmm模型聚类中的特异性指标筛选出去冗余的hmm模型

6.使用去冗余后的hmm模型与病原菌和非病原菌的蛋白数据进行比对，并根据每个模型比对上该蛋白数据的序列数量，使用tspex软件计算每个hmm模型对于病原菌和非病原菌的特异性指标，使用特异性指标作为marker筛选的阈值，符合标准的hmm模型则作为病原菌的蛋白marker和非病原的蛋白marker

7.利用每个模型比对上这些蛋白marker的数量计算出用于机器学习训练的各种特征值，利用生成的各类细菌序列中的特征矩阵使用各类机器学习模型进行分类器监督学习，得到分类
