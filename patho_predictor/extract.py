import pandas as pd
from Bio import SeqIO
import random
from tqdm import tqdm
import argparse

def extract_sequences(input_file, output_file, num_sequences):
    sequences = list(SeqIO.parse(input_file, "fasta"))

    if len(sequences) > num_sequences:
        sampled_sequences = random.sample(sequences, num_sequences)
    else:
        sampled_sequences = sequences

    with open(output_file, "w") as output_handle:
        SeqIO.write(sampled_sequences, output_handle, "fasta")

def get_files(files, num_sequences):
    for file in tqdm(files, desc="Processing Files"):
        input_file = '../data_short_read/' + file
        output_file = '../extracted_1000_short_reads/' + file
        extract_sequences(input_file, output_file, num_sequences)
        

def main():
    parser = argparse.ArgumentParser(description="Extract sequences from FASTA files")
    parser.add_argument(
        "--num-sequences", type=int, default=1000, help="Number of sequences to extract from each file"
    )
    args = parser.parse_args()

    df = pd.read_csv('data_labels.csv')

    get_files(df['filename_x'], args.num_sequences)

if __name__ == "__main__":
    main()