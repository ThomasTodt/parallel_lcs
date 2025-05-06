import random
from pathlib import Path


# Function to generate a random DNA sequence of a given length
def generate_random_dna(length):
    return "".join(random.choices("ACGT", k=length))


# Function to generate a sequence with a specific pattern repeated
def generate_sequence_with_pattern(pattern, repeats, total_length):
    pattern_seq = pattern * repeats
    if len(pattern_seq) > total_length:
        raise ValueError("Pattern sequence length exceeds total sequence length.")
    remaining_length = total_length - len(pattern_seq)
    seq = (
        generate_random_dna(remaining_length // 2)
        + pattern_seq
        + generate_random_dna(remaining_length - remaining_length // 2)
    )
    return seq


# Example usage
if __name__ == "__main__":
    pattern = "ATGC"
    sequence_size = 10

    # generate a sequences with increasing size, by 10 times each interaction
    for i in range(1, 7):
        seq1 = generate_sequence_with_pattern(pattern, i, sequence_size)

        with Path(f"seq2_{sequence_size}.txt").open("w") as f:
            f.write(seq1)

        sequence_size *= 10

        print(
            f"Generated sequence with pattern with size {sequence_size} and length {len(seq1)}",
        )
