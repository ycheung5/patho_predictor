import pandas as pd
from sklearn.model_selection import train_test_split
from joblib import Parallel, delayed
from Bio import SeqIO
import os

def data_labels_train_test_split(data_path='data_labels.csv', test_size=0.1, random_state=42, save_csv=False, train_path='train.csv', test_path='test.csv'):
    try:
        df = pd.read_csv(data_path)
    except FileNotFoundError:
        print(f"Error: The file {data_path} does not exist.")
        return None, None
    except pd.errors.EmptyDataError:
        print(f"Error: The file {data_path} is empty.")
        return None, None

    # Check if required columns are in the dataframe
    required_columns = {'filename_x', 'group'}
    if not required_columns.issubset(df.columns):
        print(f"Error: The dataframe must include the columns {required_columns}.")
        return None, None

    df_group = df[['filename_x', 'group']]
    train_df, test_df = train_test_split(df_group, test_size=test_size, stratify=df_group['group'], random_state=random_state)

    if save_csv:
        train_df.to_csv(train_path, index=False)
        test_df.to_csv(test_path, index=False)
        print(f"Train data saved to {train_path}")
        print(f"Test data saved to {test_path}")

    return train_df, test_df


def sample_by_group_fraction(kmers_df, sample_size=0.1):
    def sample_group(group):
        n_samples = max(1, int(len(group) * sample_size))  # Ensures at least one sample
        return group.sample(n=n_samples)

    # Applying the function to each group
    return kmers_df.groupby('filename_x').apply(sample_group).reset_index(drop=True)
# test_kmers_df_10 = sample_by_group_fraction(test_kmers_df)

#file_dir = "/projects/cs6824_f19/Tang_Yat_Vineeth/extracted_1000_short_reads_40/"
def process_fasta_file(filename, file_dir, group):
    file_path = os.path.join(file_dir, filename)  # Ensure file path is correctly formed
    data = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        data.append({'sequence': sequence, 'group': group})
    return pd.DataFrame(data)

def read_fasta_files(file_dir, selected_samples, n_jobs=-1):
    results = Parallel(n_jobs=n_jobs)(
        delayed(process_fasta_file)(row['filename_x'], file_dir, row['group'])
        for index, row in selected_samples.iterrows())
    return pd.concat(results, ignore_index=True)


def write_seq_group_df_to_fasta(df, filename):
    with open(filename, 'w') as file:
        for index, row in df.iterrows():
            header = f'>{index}_{row["group"]}'
            sequence = row['sequence']
            file.write(f'{header}\n{sequence}\n')

def read_seq_group_df_back_from_fasta(filename):
    sequences = []
    groups = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Extract group from the header, assuming format >index_group
                group = line.strip().split('_')[-1]
                groups.append(group)
            else:
                sequences.append(line.strip())
                
    return pd.DataFrame({'sequence': sequences, 'group': groups})

