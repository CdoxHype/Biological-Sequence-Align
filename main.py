import numpy as np
import csv
import sys

# Author: Julio Camacho Alicea
# Description: This script performs sequence alignment using dynamic programming and the Needleman-Wunsch algorithm.
# It reads two sequences from a CSV file and computes their optimal alignment, outputting the aligned sequences and their score.

def ParserFile(file):
    """
    Parses a CSV file containing two sequences to be compared, and generates a matrix for sequence alignment.

    Args:
        file (str): The path to the CSV file containing the sequences.
    """
    seq_one, seq_two = '', ''
    try:
        # Open the CSV file for reading
        with open(file, newline='') as cvsfile:
            reader = csv.reader(cvsfile)
            for row_num, row in enumerate(reader, start=1):
                # Extract sequences from the CSV file
                seq_one = row[0]
                seq_two = row[1]
                one_size = len(seq_one) + 1  # Length of first sequence including a gap
                two_size = len(seq_two) + 1  # Length of second sequence including a gap
                matrix = np.full((one_size, two_size), None)  # Initialize the matrix with None values
                
                # Populate the matrix with alignment values
                result = PopulateMatrix(one_size, two_size, seq_one, seq_two, matrix)
                print(result)  # Print the result after matrix population
                
    except Exception as e:
        print(f"An error occurred: {e}")  # Handle any errors encountered while parsing the file


def PopulateMatrix(seq_rows, seq_columns, sequence1, sequence2, matx):
    """
    Populates the matrix with sequence alignment values based on a scoring system.

    Args:
        seq_rows (int): Number of rows in the matrix (length of first sequence + 1).
        seq_columns (int): Number of columns in the matrix (length of second sequence + 1).
        sequence1 (str): The first sequence to be aligned.
        sequence2 (str): The second sequence to be aligned.
        matx (ndarray): The matrix to be populated.

    Returns:
        str: The final aligned sequences and the score.
    """
    seq1 = sequence1
    seq2 = sequence2
    matrix = matx
    gap_penalty = -2  # Gap penalty score

    # Populate the matrix with scoring values
    for row in range(seq_rows):
        for column in range(seq_columns):
            if row == 0 and column == 0:  # Base case (top-left corner)
                matrix[row, column] = 0
            elif row == 0:  # First row, fill with gap penalties
                matrix[row, column] = matrix[row, column - 1] + gap_penalty
            elif column == 0:  # First column, fill with gap penalties
                matrix[row, column] = matrix[row - 1, column] + gap_penalty
            else:
                # Calculate the three possible values (match/mismatch, deletion, insertion)
                match_mismatch = matrix[row - 1, column - 1] + ScoringMatrix(row, column, seq1, seq2)
                delete = matrix[row - 1, column] + gap_penalty
                insert = matrix[row, column - 1] + gap_penalty

                # Choose the maximum value for the matrix cell
                matrix[row, column] = max(match_mismatch, delete, insert)

    # Perform backtracking to get the optimal alignment
    return BackTracking(matrix, seq1, seq2)


def ScoringMatrix(current_row, current_column, seq1, seq2):
    """
    Computes the score for the alignment of two letters from each sequence.

    Args:
        current_row (int): The current row index.
        current_column (int): The current column index.
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        int: The score for the alignment (1 for match, -1 for mismatch).
    """
    SeqLetterOne = seq1[current_row - 1]
    SeqLetterTwo = seq2[current_column - 1]

    if SeqLetterOne == SeqLetterTwo:
        return 1  # Match: +1 score
    else:
        return -1  # Mismatch: -1 score


def BackTracking(matx, sequence1, sequence2):
    """
    Performs backtracking to generate the aligned sequences based on the populated matrix.

    Args:
        matx (ndarray): The populated scoring matrix.
        sequence1 (str): The first sequence to be aligned.
        sequence2 (str): The second sequence to be aligned.

    Returns:
        str: The final aligned sequences and the score.
    """
    seq_one = sequence1
    seq_two = sequence2
    matrix = matx

    seq_one_final = ''  # Aligned first sequence
    seq_two_final = ''  # Aligned second sequence

    # Start at the bottom-right corner of the matrix
    row_index = matrix.shape[0] - 1
    column_index = matrix.shape[1] - 1

    # Perform backtracking to reconstruct the alignment
    while row_index > 0 and column_index > 0:
        # Calculate the possible scores for diagonal, up, and left moves
        diagonal_score = int(matrix[row_index - 1, column_index - 1] + ScoringMatrix(row_index, column_index, seq_one, seq_two))
        upper_score = int(matrix[row_index - 1, column_index] + -2)  # Gap penalty
        side_score = int(matrix[row_index, column_index - 1] + -2)  # Gap penalty

        # Determine the direction to move (diagonal, up, or left)
        if diagonal_score > upper_score and diagonal_score > side_score:
            seq_one_final = seq_one[row_index - 1] + seq_one_final
            seq_two_final = seq_two[column_index - 1] + seq_two_final
            row_index -= 1
            column_index -= 1
        elif upper_score > side_score:
            seq_one_final = seq_one[row_index - 1] + seq_one_final
            seq_two_final = '-' + seq_two_final  # Insert gap in sequence 2
            row_index -= 1
        else:
            seq_one_final = '-' + seq_one_final  # Insert gap in sequence 1
            seq_two_final = seq_two[column_index - 1] + seq_two_final
            column_index -= 1

    # Handle cases where one sequence has been fully aligned
    if column_index == 0:
        while row_index != 0:
            seq_one_final = seq_one[row_index - 1] + seq_one_final
            seq_two_final = '-' + seq_two_final
            row_index -= 1
    elif row_index == 0:
        while column_index != 0:
            seq_two_final = seq_two[column_index - 1] + seq_two_final
            seq_one_final = '-' + seq_one_final
            column_index -= 1

    # Output the final aligned sequences and the score
    print(f'{seq_one_final} {seq_two_final} {matrix[-1, -1]}')
    return f'{seq_one_final} {seq_two_final} {matrix[-1, -1]}'


def main():
    """
    Main function that initiates the sequence alignment process by calling the ParserFile function.

    Reads the sequences from a CSV file and prints the result of the alignment.
    """
    ParserFile('./sequences.csv')  # Path to the CSV file containing sequences

if __name__ == "__main__":
    main()  # Execute the main function when the script is run
