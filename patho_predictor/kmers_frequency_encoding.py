import argparse
import pandas as pd
from Bio import SeqIO
import os
from itertools import product
from joblib import Parallel, delayed
from tqdm import tqdm

def encode_sequence(sequence, k):
    all_kmers = [''.join(p) for p in product('ATCG', repeat=k)]
    kmer_counts = {kmer: 0 for kmer in all_kmers}
    sequence = ''.join(char for char in sequence if char in ['A', 'T', 'C', 'G'])
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
    encoded_df = pd.DataFrame([kmer_counts.values()], columns=all_kmers)
    sum_counts = encoded_df.sum(axis=1)
    encoded_df = encoded_df.divide(sum_counts, axis=0)
    encoded_df = encoded_df.fillna(0)
    return encoded_df

# k_mers = [3]
# file_dir = '../data'
def read_fasta_to_kmers(df, file_dir, k_mers):
    def process_file(filename, group):
        dfs_per_sequence = []
        file_path = os.path.join(file_dir, filename)
        for record in SeqIO.parse(file_path, "fasta"):
            sequence = str(record.seq)
            encoded_sequence_dfs = [encode_sequence(sequence, k) for k in k_mers]
            merged_sequence_df = pd.concat(encoded_sequence_dfs, axis=1)
            merged_sequence_df['filename_x'] = filename
            merged_sequence_df['group'] = group
            dfs_per_sequence.append(merged_sequence_df)
        return dfs_per_sequence
    results = Parallel(n_jobs=-1)(delayed(process_file)(filename, group) for filename, group in tqdm(zip(df['filename_x'], df['group']), total=df.shape[0]))
    final_df = pd.concat([pd.concat(file_dfs, ignore_index=True) for file_dfs in results], ignore_index=True)
    final_df = final_df.fillna(0)
    return final_df

def main(args):
    df = pd.read_csv(args.input_data_labels)
    for k in args.kmer_lengths:
        print(f"Processing k-mer length: {k}")
        kmers_df = read_fasta_to_kmers(df, args.file_dir, [k])
        output_path = f"{args.output_name}_{k}mers.csv"
        kmers_df.to_csv(output_path, index=False)
        print(f"Encoded k-mer data saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Encode k-mer features from fasta files.")
    parser.add_argument("--kmer_lengths", type=int, nargs='+', required=True, help="List of k-mer lengths to analyze")
    parser.add_argument("--input_data_labels", type=str, default="data_labels.csv", help="CSV file with data labels")
    parser.add_argument("--file_dir", type=str, default="../data", help="Directory with fasta files")
    parser.add_argument("--output_name", type=str, default="data_labels", help="Base name for output CSV files")
    args = parser.parse_args()

    main(args)
# nohup python3 kmers_frequency_encoding.py --kmer_lengths "3" --input_data_labels "train.csv" --file_dir "/projects/cs6824_f
# 19/Tang_Yat_Vineeth/extracted_1000_short_reads_40/" --output_name "train_3mers_extracted_sr_40" > 3mers_extract_40.log &
# nohup python3 kmers_frequency_encoding.py --kmer_lengths "3" --input_data_labels "train.csv" --file_dir "/projects/cs6824_f
# 19/Tang_Yat_Vineeth/extracted_1000_short_reads/" --output_name "train_3mers_extracted_sr" > 3mers_extract.log &