import pandas as pd
import sys

def calculate_percentage_of_ones(df):
    results = []
    for col in df.columns:
        if col.endswith('-g'):
            total_ones = df[col].sum()
            total_entries = len(df[col])
            percentage_ones = total_ones / total_entries
            results.append((col, percentage_ones))
    return results

def main(input_file):
    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # Calculate percentage of 1s for columns ending with '-g'
    results = calculate_percentage_of_ones(df)
    
    # Output results in TSV format
    for col, percentage in results:
        print(f"{input_file.split('.')[0]}\t{col}\t{percentage:.2f}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input_file.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    main(input_file)

