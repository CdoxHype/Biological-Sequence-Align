# Needleman-Wunsch Algorithm

This project implements the **Needleman-Wunsch** algorithm, which is used for sequence alignment. The algorithm is widely used in bioinformatics for comparing biological sequences, such as DNA, RNA, or protein sequences, to find the best possible alignment.

## Description

The Needleman-Wunsch algorithm is a dynamic programming algorithm that aligns two sequences by optimizing a score based on a scoring matrix. This project includes:

- **Parser**: Reads sequences from a CSV file.
- **Matrix Population**: Fills a matrix with values based on the sequences.
- **Backtracking**: Finds the optimal alignment by tracing back through the matrix.
- **Scoring**: Implements a scoring system for sequence matching and gap penalties.

## Features

- Reads sequence data from a CSV file.
- Performs global sequence alignment using the Needleman-Wunsch algorithm.
- Displays the aligned sequences and the final score.
- Customizable scoring matrix and gap penalties.

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/your-username/needleman-wunsch.git
    ```

2. Navigate into the project directory:
    ```bash
    cd needleman-wunsch
    ```

## Usage

1. Make sure you have a CSV file (`sequences.csv`) with two sequences to align. The file should have two columns with each row containing one pair of sequences. Example:

    ```csv
    AGCTG,AGCT
    ACGT,ACGT
    ```

2. Run the script:
    ```bash
    python needleman_wunsch.py
    ```

   This will output the aligned sequences and the alignment score.

## Functions

- **ParserFile(file)**: Reads a CSV file and extracts two sequences.
- **PopulateMatrix(seq_rows, seq_columns, sequence1, sequence2, matx)**: Populates the alignment matrix based on the sequences.
- **ScoringMatrix(current_row, current_column, seq1, seq2)**: Computes the score for sequence matching.
- **BackTracking(matx, sequence1, sequence2)**: Traces back through the matrix to find the optimal alignment.

## Author

**Julio Camacho Alicea**
