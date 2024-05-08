# Pathogenicity of DNA Sequences Predictor

This tool predicts the pathogenicity of sequences read from a FASTA file using a pre-trained RandomForest model. The sequences are encoded into k-mer frequencies, and the model outputs the probability of each sequence being pathogenic.

## Features
 
- Reads sequences from a FASTA file.
- Encodes sequences into k-mer frequencies.
- Uses a RandomForest model to predict pathogenicity.
- Outputs predictions to a CSV file.

## Requirements

- Python 3.6+
- pandas
- Biopython
- joblib

## Installation

First, clone the repository or download the source code. Then, navigate to the project directory and install the required Python packages using pip:

```bash
pip install -r requirements.txt
```

## Usage

To run the tool, you'll need to specify the path to the pre-trained RandomForest model, the input FASTA file, and the desired output CSV file for the predictions.

```bash
python patho_predict.py --model_path "path/to/model.joblib" --input_fasta "path/to/input.fasta" --output_csv "path/to/output.csv"
```

Parameters
--model_path (required): Path to the pre-trained RandomForest model file (.joblib).
--input_fasta (required): Path to the FASTA file containing sequences to predict.
--output_csv (required): Path where the prediction output CSV file will be saved.

## Output
The output CSV file will contain three columns:

sequence_id: The ID of the sequence from the FASTA file.
prediction_value: The probability of the sequence being pathogenic.
label: A label indicating "pathogenic" if the probability is 0.5 or higher, and "non-pathogenic" otherwise.

## Demo

Run the script with the following command:

```bash
python patho_predict.py --model_path "demo/random_forest.joblib" --input_fasta "demo/demo.fna" --output_csv "result_demo.csv"
```
