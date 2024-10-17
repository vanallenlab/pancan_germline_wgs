import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def plot_qc_metrics(pre_filtered_df, post_filtered_df):
    # Check if all columns of interest are present in both DataFrames
    columns_of_interest = ['n_grafpop_snps', 'hq_het_rate', 'charr', 'inconsistent_ab_het_rate',
                           'median_coverage', 'chrX_ploidy', 'chrY_ploidy']

    missing_columns = []
    for col in columns_of_interest:
        if col not in pre_filtered_df.columns:
            missing_columns.append((args.pre_filtered, col))
        if col not in post_filtered_df.columns:
            missing_columns.append((args.post_filtered, col))

    if missing_columns:
        for filename, column in missing_columns:
            print(f"Error: {filename} doesn't contain the '{column}' column.")
        return

    # Filter 'Post-filtered' DataFrame for rows where global_qc_pass is TRUE
    print(post_filtered_df.head())
    post_filtered_df = post_filtered_df[post_filtered_df['global_qc_pass'] == True]
    print(post_filtered_df[post_filtered_df['global_qc_pass'] == True].head())

    # Define columns of interest
    for metric in columns_of_interest:
        # Create a violin plot for pre-filtered and post-filtered datasets
        plt.figure(figsize=(10, 6))

        # Violin plot
        plt.subplot(1, 2, 1)  # Subplot for violin plot
        sns.violinplot(data=[pre_filtered_df[metric], post_filtered_df[metric]], palette="Set3")
        plt.title(f'{metric} - Quality Control')
        plt.xlabel('Dataset')
        plt.ylabel(metric)
        plt.xticks([0, 1], ['Pre-filtered', 'Post-filtered'])

        # Box plot
        plt.subplot(1, 2, 2)  # Subplot for box plot
        sns.boxplot(data=[pre_filtered_df[metric], post_filtered_df[metric]], palette="Set3")
        plt.title(f'{metric} - Quality Control')
        plt.xlabel('Dataset')
        plt.ylabel(metric)
        plt.xticks([0, 1], ['Pre-filtered', 'Post-filtered'])

        # Adjust layout and spacing
        plt.tight_layout()

        # Save the plot
        plt.savefig(f'{metric}_qc.png')
        plt.close()


def parse_arguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Plot quality control metrics")
    # Add arguments
    parser.add_argument("--pre-filtered", required=True, help="Path to pre-filtered TSV file")
    parser.add_argument("--post-filtered", required=True, help="Path to post-filtered TSV file")
    # Parse arguments
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Read pre-filtered and post-filtered TSV files into DataFrames
    try:
        pre_filtered_df = pd.read_csv(args.pre_filtered, sep='\t')
        post_filtered_df = pd.read_csv(args.post_filtered, sep='\t')
    except FileNotFoundError as e:
        print(f"Error: {e.filename} not found.")
        exit(1)

    # Call plot function with DataFrame arguments
    plot_qc_metrics(pre_filtered_df, post_filtered_df)