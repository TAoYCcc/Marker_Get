class Config:
    def __init__(self):
        self.input_path = r'/home/huhaoyu/NHP_ftp.csv'
        self.output_path = r'/home/huhaoyu/NHP/'
        self.cpu_worker_num = 180
        # file suffix of final url
        self.file_type1 = '_genomic.fna.gz'
        self.file_type2 = '_protein.faa.gz'
        self.csv_col = 'ftp_path'

        # times of re-downloads
        self.retry_times = 5
        self.fail_log_path = r'/home/huhaoyu/NHP/logNHP.txt'
