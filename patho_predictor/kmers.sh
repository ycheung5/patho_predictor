#!/bin/bash
#SBATCH --cpus-per-task=100
#SBATCH --time=100:00:00
#SBATCH --partition=normal_q
#SBATCH --account=cs3824_fall23

nohup python3 kmers_frequency_encoding.py --kmer_lengths "3" --input_data_labels "train.csv" --file_dir "/projects/cs6824_f19/Tang_Yat_Vineeth/extracted_1000_short_reads_40/" --output_name "train_3mers_extracted_sr_40" > 3mers_extract_40_.log 2>&1 & 

nohup python3 kmers_frequency_encoding.py --kmer_lengths "3" --input_data_labels "test.csv" --file_dir "/projects/cs6824_f19/Tang_Yat_Vineeth/extracted_1000_short_reads_40/" --output_name "test_3mers_extracted_sr_40" > 3mers_extract_40_test_.log 2>&1 & 

nohup python3 kmers_frequency_encoding.py --kmer_lengths "3" --input_data_labels "train.csv" --file_dir "/projects/cs6824_f19/Tang_Yat_Vineeth/extracted_1000_short_reads/" --output_name "train_3mers_extracted_sr" > 3mers_extract_.log 2>&1 & 
