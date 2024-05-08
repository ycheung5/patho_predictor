import pandas as pd
import numpy as np
from collections import Counter
import math
from itertools import combinations, product


# Define the feature extraction functions
def clean_sequence(sequence):
    return ''.join([nuc for nuc in sequence if nuc in 'ATCG'])

def gc_content(sequence):
    sequence = clean_sequence(sequence)
    return (sequence.count('G') + sequence.count('C')) / len(sequence) if sequence else 0

def shannon_entropy(sequence):
    sequence = clean_sequence(sequence)
    count = Counter(sequence)
    total = len(sequence)
    return -sum((freq / total) * math.log2(freq / total) for freq in count.values()) if total > 0 else 0

def nucleotide_frequency(sequence):
    sequence = clean_sequence(sequence)
    freq = Counter(sequence)
    total = len(sequence) if sequence else 1
    return [freq.get(nuc, 0) / total for nuc in 'ATCG']

def binary_encoding(sequence):
    encoding_map = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1]}
    sequence = clean_sequence(sequence)
    encoded = [encoding_map[nuc] for nuc in sequence]
    return np.array(encoded).flatten().tolist()[:100]  # Flatten and truncate for simplicity

def fourier_transform(sequence):
    sequence = clean_sequence(sequence)
    numeric_seq = np.array([{'A': 0, 'C': 1, 'G': 2, 'T': 3}.get(nuc, 0) for nuc in sequence])
    transformed = np.fft.fft(numeric_seq)
    return np.abs(transformed)[:10].tolist()  # Return first 10 elements of the Fourier transform

def kmer_frequencies(sequence, k=3, alphabet='ACGT'):
    sequence = clean_sequence(sequence)
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    kmer_counts = Counter(kmers)
    total_kmers = sum(kmer_counts.values())
    
    # Generate all possible k-mers from the given alphabet, sorted to ensure consistent order
    all_kmers = [''.join(p) for p in product(alphabet, repeat=k)]
    all_kmers.sort()  # Ensure the order is consistent, although product should already do this in lexicographical order

    # Initialize the frequency dictionary with all possible k-mers set to 0
    kmer_freq = {kmer: 0 for kmer in all_kmers}
    
    # Update frequencies for k-mers found in the sequence
    for kmer, count in kmer_counts.items():
        if kmer in kmer_freq:  # Check if kmer is valid per the given alphabet and length
            kmer_freq[kmer] = count / total_kmers

    return kmer_freq



def dinucleotide_properties(sequence):
    sequence = clean_sequence(sequence)
    properties = {}
    for i in range(len(sequence) - 1):
        dinucleotide = sequence[i:i+2]
        properties[dinucleotide] = properties.get(dinucleotide, 0) + 1
    return properties

def positional_nucleotide_frequencies(sequence):
    sequence = clean_sequence(sequence)
    pos_freq = {f"pos_{i}_{nuc}": 0 for i in range(10) for nuc in 'ATCG'}  # Example for first 10 positions
    for i, nuc in enumerate(sequence[:10]):  # Only first 10 positions considered for simplicity
        pos_freq[f"pos_{i}_{nuc}"] += 1
    return pos_freq

def kmer_type_1(sequence, k=2):
    # Normalized k-mer frequencies
    sequence = clean_sequence(sequence)
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)
    total_kmers = sum(kmer_counts.values())
    return {kmer: count / total_kmers for kmer, count in kmer_counts.items()}

def mismatch_profile(sequence, k=2, mismatch=1):
    # This requires a more complex implementation that includes allowing for mismatches
    # Simplified version: count kmers with up to 'mismatch' mismatches
    from itertools import product
    sequence = clean_sequence(sequence)
    kmers = set([sequence[i:i+k] for i in range(len(sequence)-k+1)])
    mismatch_kmers = set()
    
    def generate_mismatches(kmer):
        bases = ['A', 'C', 'G', 'T']
        for positions in combinations(range(k), mismatch):
            for replacements in product(bases, repeat=mismatch):
                new_kmer = list(kmer)
                for pos, rep in zip(positions, replacements):
                    new_kmer[pos] = rep
                mismatch_kmers.add(''.join(new_kmer))
    
    for kmer in kmers:
        generate_mismatches(kmer)
    
    # Filter to keep only those that are within the sequence
    valid_kmers = [mk for mk in mismatch_kmers if mk in kmers]
    return len(valid_kmers) / len(kmers) if kmers else 0

def nucleotide_auto_covariance(sequence, lag=1):
    sequence = clean_sequence(sequence)
    numeric_seq = [1 if x == 'A' else 2 if x == 'C' else 3 if x == 'G' else 4 for x in sequence]
    if len(numeric_seq) <= lag:
        return 0
    mean = sum(numeric_seq) / len(numeric_seq)
    covariances = [(numeric_seq[i] - mean) * (numeric_seq[i+lag] - mean) for i in range(len(numeric_seq) - lag)]
    return sum(covariances) / len(covariances)

# def z_curve(sequence):
#     sequence = clean_sequence(sequence)
#     x = y = z = 0
#     x_values, y_values, z_values = [], [], []

#     for nuc in sequence:
#         if nuc == 'A':
#             x += 1
#         elif nuc == 'C':
#             y += 1
#         elif nuc == 'G':
#             z += 1
#         elif nuc == 'T':
#             x -= 1
#             y -= 1
#             z -= 1
#         x_values.append(x)
#         y_values.append(y)
#         z_values.append(z)

#     return x_values, y_values, z_values


def pssm(sequence, motif_length=3):
    from itertools import product
    sequence = clean_sequence(sequence)
    kmers = [''.join(k) for k in product('ATCG', repeat=motif_length)]
    pssm_scores = {k: 0 for k in kmers}
    total_kmers = len(sequence) - motif_length + 1

    for i in range(total_kmers):
        kmer = sequence[i:i + motif_length]
        if kmer in pssm_scores:
            pssm_scores[kmer] += 1

    for kmer in pssm_scores:
        pssm_scores[kmer] /= total_kmers

    return pssm_scores

# def accumulated_nucleotide_frequency(sequence):
#     sequence = clean_sequence(sequence)
#     freq = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
#     accumulated_freqs = []

#     total = 0
#     for nuc in sequence:
#         freq[nuc] += 1
#         total += 1
#         accumulated_freqs.append({nuc: freq[nuc] / total for nuc in 'ACGT'})

#     return accumulated_freqs

features_dict = {
    'GC_Content': gc_content,
    'Shannon_Entropy': shannon_entropy,
    'Nucleotide_Frequency': nucleotide_frequency,
    'Binary_Encoding': binary_encoding,
    'Fourier_Transform': fourier_transform,
    'Kmer_Frequencies': kmer_frequencies,
    'Dinucleotide_Properties': dinucleotide_properties,
    'Positional_Nucleotide_Frequencies': positional_nucleotide_frequencies,
    'Kmer_Type_1': lambda seq: kmer_type_1(seq, k=2),
    'Mismatch_Profile': lambda seq: mismatch_profile(seq, k=2, mismatch=1),
    'Nucleotide_Auto_Covariance': lambda seq: nucleotide_auto_covariance(seq, lag=1),
    #'Z_Curve': z_curve,
    'PSSM': lambda seq: pssm(seq, motif_length=3),
    #'Accumulated_Nucleotide_Frequency': accumulated_nucleotide_frequency
}
    
# Function to apply selected features and expand list outputs
def apply_selected_features(df, feature_names):
    for feature in feature_names:
        if feature in features_dict:
            feature_data = df['sequence'].apply(features_dict[feature])
            if isinstance(feature_data.iloc[0], dict):  # Check if the output is a dictionary
                feature_df = pd.json_normalize(feature_data)
                feature_df.columns = [f"{feature}_{col}" for col in feature_df.columns]
                df = pd.concat([df, feature_df], axis=1)
            elif isinstance(feature_data.iloc[0], list):  # Check if the output is a list
                feature_df = pd.DataFrame(feature_data.tolist(), columns=[f"{feature}_{i}" for i in range(len(feature_data.iloc[0]))])
                df = pd.concat([df, feature_df], axis=1)
            else:
                df[feature] = feature_data
        else:
            print(f"Feature {feature} not recognized.")
    return df

# selected_features = ['GC_Content']
# df = apply_selected_features(train_df, selected_features)

def process_fasta(fasta_file, selected_features):
    # Process the fasta file and extract features
    sequences = [record.seq for record in SeqIO.parse(fasta_file, 'fasta')]
    df = pd.DataFrame({'sequence': sequences})

    # Apply selected features
    df = apply_selected_features(df, selected_features)
    return df


def user_interaction():
    # User interaction to select features and process files
    fasta_file = input("Please enter the path to your FASTA file: ")
    print("Available features: ", list(features_dict.keys()))
    selected_features = input("Enter the features you want to apply (comma separated): ").split(',')
    selected_features = [feature.strip() for feature in selected_features if feature.strip() in features_dict]

    df = process_fasta(fasta_file, selected_features)
    output_file = input("Enter the name for the output CSV file: ")
    df.to_csv(output_file, index=False)
    print(f"Data has been encoded and saved to {output_file}")

if __name__ == "__main__":
    user_interaction()