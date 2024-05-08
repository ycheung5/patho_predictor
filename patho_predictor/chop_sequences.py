import os
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed 


def chop_single_file(input_file, output_dir, chunk_size=150):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create the output file path
    filename = os.path.basename(input_file)
    output_path = os.path.join(output_dir, filename)
    
    # Open the input file and create the output file
    with open(input_file, "r") as input_handle, open(output_path, "w") as output_handle:
        # Parse the sequences in the FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence_length = len(record.seq)
            full_chunks = sequence_length // chunk_size
            
            # Chop the sequence into full chunks and write each to the output file
            for i in range(full_chunks):
                start = i * chunk_size
                end = start + chunk_size
                new_record = record[start:end]
                new_record.id = f"{record.id}_chunk{i+1}"
                new_record.description = f"chunk {i+1} of {record.id}"
                SeqIO.write(new_record, output_handle, "fasta")
                
# input_file_path = '../data/GCF_000005825.2_ASM582v2_genomic.fna'
# output_directory = '../data_short_read'
# chop_single_file(input_file_path, output_directory)

def process_file(filename, input_dir, output_dir, overlap, chunk_size=150):
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, filename)
    
    # Open the input file and create the output file
    with open(input_path, "r") as input_handle, open(output_path, "w") as output_handle:
        # Parse the sequences in the FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence_length = len(record.seq)
#             full_chunks = sequence_length // chunk_size

            # Chop the sequence into full chunks and write each to the output file
#             for i in range(full_chunks):
#                 start = i * chunk_size
#                 end = start + chunk_size
#                 new_record = record[start:end]
#                 new_record.id = f"{record.id}_chunk{i+1}"
#                 new_record.description = f"chunk {i+1} of {record.id}"
#                 SeqIO.write(new_record, output_handle, "fasta")
            curr = 0
            i = 0
            while curr < sequence_length:
                start = curr
                end = curr + chunk_size
                if sequence_length < end:
                    break
                new_record = record[start:end]
                new_record.id = f"{record.id}_chunk{i+1}"
                new_record.description = f"chunk {i+1} of {record.id}"
                SeqIO.write(new_record, output_handle, "fasta")
                curr += (chunk_size - overlap)
                i+=1

def chop_sequences_parallel(input_dir, output_dir, overlap=0, chunk_size=150):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Get all fasta files in the input directory
    files = [f for f in os.listdir(input_dir) if f.endswith('.fna')]
    
    # Set up a process pool to execute processes in parallel
    with ProcessPoolExecutor() as executor:
        # Map the process_file function to the files
        futures = [executor.submit(process_file, f, input_dir, output_dir, overlap, chunk_size) for f in files]
        
        # Wait for all futures to complete
        for future in futures:
            future.result()  # Getting the result to handle exceptions

def count_sequences_in_file(file_path):
    """ Helper function to count sequences in a single file """
    count = 0
    with open(file_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            count += 1
    return count

def count_sequences_in_fasta(directory):
    # Initialize a list to store the results
    results = []

    # Prepare a list of tasks with file paths
    file_paths = [os.path.join(directory, filename)
                  for filename in os.listdir(directory)
                  if filename.endswith(".fasta") or filename.endswith(".fna")]

    # Use ProcessPoolExecutor to handle parallel processing
    with ProcessPoolExecutor() as executor:
        future_to_file = {executor.submit(count_sequences_in_file, file_path): file_path
                          for file_path in file_paths}

        for future in as_completed(future_to_file):  # Corrected reference here
            file_path = future_to_file[future]
            filename = os.path.basename(file_path)
            count = future.result()
            results.append({'Filename': filename, 'Sequence Count': count})

    # Create a DataFrame from the collected data
    df = pd.DataFrame(results)
    return df

def main():
    parser = argparse.ArgumentParser(description='Chop sequences into smaller chunks with overlap.')
    parser.add_argument('-i', '--input_dir', type=str, default='../data/', help='Input directory containing fasta files.')
    parser.add_argument('-o', '--output_dir', type=str, default='../data_short_read_40/', help='Output directory for chunked fasta files.')
    parser.add_argument('-l', '--overlap', type=int, default=40, help='Overlap between chunks.')
    parser.add_argument('-c', '--chunk_size', type=int, default=150, help='Size of each chunk.')

    args = parser.parse_args()

    chop_sequences_parallel(args.input_dir, args.output_dir, args.overlap, args.chunk_size)

if __name__ == "__main__":
    main()