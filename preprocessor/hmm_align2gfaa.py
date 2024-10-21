import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


def hmmsearch_task(hmm_file, sequence_db, res_folder, base_name, db_name):
    output_file = os.path.join(res_folder, base_name, f"{base_name}_to_{db_name}_results.txt")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    cmd = ['hmmsearch', '--cpu', '10', '-E', '1e-5', hmm_file, sequence_db]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
        stdout, stderr = proc.communicate()
        if proc.returncode == 0:
            with open(output_file, 'w') as f:
                f.write(stdout)
        else:
            print(f"Error running hmmsearch for {hmm_file} and {sequence_db}: {stderr}")


def main():
    hmm_folder = "/home/taoyechen/Balance_process/hmm_clustered"
    sequence_db_folder = "/home/taoyechen/NCBI_oneS_balance/HP_faa_296"
    res_folder = "/home/taoyechen/Balance_process/HP2hmmC"
    concurrency = 20

    hmm_files = [os.path.join(hmm_folder, f) for f in os.listdir(hmm_folder) if f.endswith('.hmm')]
    futures = []
    with ThreadPoolExecutor(max_workers=concurrency) as executor:
        for hmm_file in hmm_files:
            base_name = os.path.splitext(os.path.basename(hmm_file))[0]
            sequence_dbs = [os.path.join(sequence_db_folder, f) for f in os.listdir(sequence_db_folder) if
                            f.endswith('.faa')]
            for sequence_db in sequence_dbs:
                db_name = os.path.splitext(os.path.basename(sequence_db))[0]
                future = executor.submit(hmmsearch_task, hmm_file, sequence_db, res_folder, base_name, db_name)
                futures.append(future)

                # 等待直到达到并发数或所有任务都已提交
                if len(futures) >= concurrency:
                    for future in as_completed(futures):
                        futures.remove(future)

                        # 等待剩余的任务完成（如果有的话）
    for future in as_completed(futures):
        future.result()  # 这将引发异常（如果有的话），但在这里我们可能不想捕获它们


if __name__ == "__main__":
    main()

# 注意：上面的脚本在达到并发数时不会完全阻塞新任务的提交，
# 但它会等待一些已提交的任务完成，以便为新任务腾出空间。
# 如果需要严格的并发控制，你可能需要实现一个更复杂的任务队列系统。

# 另外，请注意，`future.result()`在这里主要是为了演示，
# 它实际上会阻塞直到任务完成。但在上面的循环中，
# 我们已经通过`as_completed`处理了已完成的任务，
# 所以这里调用`result()`可能不是必需的，除非你需要检查异常。