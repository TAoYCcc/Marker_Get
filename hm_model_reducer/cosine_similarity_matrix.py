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