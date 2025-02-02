import argparse
import pandas as pd
from Bio import SeqIO
import os
from itertools import product
from joblib import load

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

def read_fasta_to_kmers(file_path, k_mers):
    dfs_per_sequence = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        encoded_sequence_dfs = [encode_sequence(sequence, k) for k in k_mers]
        merged_sequence_df = pd.concat(encoded_sequence_dfs, axis=1)
        merged_sequence_df.insert(0, 'sequence_id', record.id)
        dfs_per_sequence.append(merged_sequence_df)
    final_df = pd.concat(dfs_per_sequence, ignore_index=True) if dfs_per_sequence else pd.DataFrame()
    return final_df

def patho_predict(model_path, input_fasta, output_csv, k_mers=[3]):
    df = read_fasta_to_kmers(input_fasta, k_mers)
    model = load(model_path)
    y_pred = model.predict_proba(df.iloc[:, 1:])[:, 1]  # Predict probabilities for the positive class
    results_df = pd.DataFrame({
        'sequence_id': df['sequence_id'],
        'prediction_value': y_pred,
        'label': ['pathogenic' if pred >= 0.5 else 'non-pathogenic' for pred in y_pred]
    })
    results_df.to_csv(output_csv, index=False)

def main():
    parser = argparse.ArgumentParser(description="Predict pathogenicity of sequences using a trained RandomForest model.")
    parser.add_argument('--model_path', type=str, required=True, help='Path to the trained model.')
    parser.add_argument('--input_fasta', type=str, required=True, help='Path to the input FASTA file.')
    parser.add_argument('--output_csv', type=str, required=True, help='Path to save the output CSV file.')
    args = parser.parse_args()
    
    patho_predict(args.model_path, args.input_fasta, args.output_csv)

if __name__ == "__main__":
    main()
