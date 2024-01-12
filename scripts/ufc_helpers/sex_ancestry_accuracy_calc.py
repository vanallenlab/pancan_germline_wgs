import pandas as pd

def process_data(df):
    # Task 1: Calculate the accuracy prediction for chrX_count, chrY_count, and SEX
    condition1 = (df['chrX_count'] == 1) & (df['chrY_count'] == 1) & (df['Sex'] == 'MALE')
    condition2 = (df['chrX_count'] == 2) & (df['chrY_count'] == 0) & (df['Sex'] == 'FEMALE')
    
    # Omit rows not meeting the conditions
    filtered_df = df[condition1 | condition2]
    
    correct_rows = filtered_df.shape[0]
    incorrect_rows = df.shape[0] - correct_rows
    
    accuracy = correct_rows / (correct_rows + incorrect_rows)
    
    print(f"\nAccuracy Prediction for chrX_count, chrY_count, and SEX: {accuracy:.2%}")
    print(f"Instances omitted: {incorrect_rows}")

    # Task 2: Calculate accuracy for all Population and Ancestry pairings
    df['Ancestry'] = df['Population'].astype(str) + '_' + df['ancestry'].astype(str)
    pairings_count = df.groupby(['ancestry', 'Population']).size().reset_index(name='Count')

    print("\nCount of Ancestry Population pairings:")
    print(pairings_count)

    return accuracy

def main():
    # Replace 'your_data.tsv' with the actual file path
    file_path = '/Users/noah/Desktop/test_sample.tsv'

    # Read the TSV file into a DataFrame
    df = pd.read_csv(file_path, sep='\t')

    # Task for each individual cohort
    cohorts = df['Cohort'].unique()
    for cohort in cohorts:
        print(f"\nProcessing cohort: {cohort}")
        cohort_df = df[df['Cohort'] == cohort]
        process_data(cohort_df)

    # Task for the combined cohort
    print("\nProcessing combined cohort:")
    process_data(df)

if __name__ == "__main__":
    main()
