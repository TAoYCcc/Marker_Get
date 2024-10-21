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