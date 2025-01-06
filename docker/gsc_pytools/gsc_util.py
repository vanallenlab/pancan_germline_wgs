import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact
import numpy as np
from firthlogist import FirthLogisticRegression
import scipy.stats as stats


# Assuming FirthLogisticRegression is already imported or defined elsewhere
# You may need to implement this if it's not available.

def logistic_regression_with_fallback(df, cancer_type, germline_event, somatic_gene, covariates=['male', 'pca_1', 'pca_2', 'pca_3', 'pca_4']):
    print(f"Cancer Type: {cancer_type}\n Germline Event: {germline_event} \n Somatic Gene: {somatic_gene}")
    
    # Filter based on cancer_type if it's not Pancancer
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

    """
    This next chunk of code is to prevent convergence failure in pancancer, when the germline event
    is not present in a cancer type. This throws off the covariates for pancancer appraoch.
    We need to keep at least 3 cancer types otherwise we throw them all out for convergence purposes.
    """
    if cancer_type == "Pancancer":
        columns_to_check = ['Breast_Diagnosis', 'Colorectal_Diagnosis', 'Kidney_Diagnosis', 'Lung_Diagnosis', 'Prostate_Diagnosis']
        
        # Identify columns to keep (those with sum != 0)
        valid_columns = [col for col in columns_to_check if df[col].sum() != 0]
        
        # Check if at least 3 of the 5 columns remain
        if len(valid_columns) < 3:
            # Drop all specified columns if fewer than 3 are valid
            df = df.drop(columns=columns_to_check)
        else:
            # Keep only the valid columns
            df = df[valid_columns + [col for col in df.columns if col not in columns_to_check]]


    # Check if the germline_event and somatic_gene are in the dataframe columns
    if germline_event not in df.columns or somatic_gene not in df.columns:
        print(f"Combination {germline_event} - {somatic_gene} not found in DataFrame")
        return
    
    # If either column has no positive cases, skip
    if df[germline_event].sum() == 0 or df[somatic_gene].sum() == 0:
        return

    # Drop rows with NaN values in the columns of interest
    df = df.dropna(subset=[germline_event, somatic_gene] + covariates)




    # Try regular logistic regression
    try:
        # Prepare the predictor (X) and response (y) variables
        X = df[[germline_event] + covariates]
        y = df[somatic_gene]

        # Normalize somatic_gene values: treat values > 1 as 1
        y = y.apply(lambda x: 1 if x > 1 else x)

        # Add a constant to the predictors for logistic regression
        X = sm.add_constant(X)

        model = sm.Logit(y, X)
        result = model.fit(disp=False)
        
        # If the model converged, extract p-value, odds ratio, and confidence intervals
        if result.mle_retvals['converged']:
            print("Regular logistic regression converged successfully!")
            
            # Extract p-value and odds ratio for the germline_gene predictor
            p_value = result.pvalues[germline_event]  # Assuming germline_event is the first predictor after the constant
            odds_ratio = np.exp(result.params[germline_event])

            # Calculate the 95% confidence interval for the odds ratio
            coef = result.params[germline_event]
            std_err = result.bse[germline_event]

            # The confidence interval for the coefficient
            conf_int_coef = [coef - 1.96 * std_err, coef + 1.96 * std_err]

            # Convert the confidence interval from log odds to odds ratio
            conf_int_or = np.exp(conf_int_coef)

            # Return the results: odds ratio, p-value, and confidence intervals
            return odds_ratio, p_value, conf_int_or[0], conf_int_or[1]
        
        else:
            raise ValueError("Regular logistic regression failed to converge.")
    
    # Fallback to Firth logistic regression if regular logistic regression fails
    except Exception as e:
        print(f"Regular logistic regression failed: {e}")
        print("Falling back to Firth logistic regression...")
        
        try:
            # Prepare the predictor (X) and response (y) variables
            X = df[[germline_event] + covariates]
            y = df[somatic_gene]

            # Normalize somatic_gene values: treat values > 1 as 1
            y = y.apply(lambda x: 1 if x > 1 else x)

            # Fit Firth logistic regression (Assuming FirthLogisticRegression is defined)
            model = FirthLogisticRegression()
            model.fit(X, y)
            
            # Extract p-value and odds ratio for the germline_gene predictor
            p_value = model.pvals_[0]  # Assuming the predictor is the second column after the constant
            odds_ratio = np.exp(model.coef_[0])
            
            # Calculate confidence intervals
            coef = model.coef_[0]
            std_err = model.bse_[0]
            conf_int_coef = [coef - 1.96 * std_err, coef + 1.96 * std_err]
            conf_int_or = np.exp(conf_int_coef)

            return odds_ratio, p_value, conf_int_or[0], conf_int_or[1]

        except Exception as e:
            print(f"Firth logistic regression failed: {e}")
            return


def firth_logistic_regression(df, cancer_type, germline_event, somatic_gene,covariates=['male','pca_1','pca_2','pca_3','pca_4']):
    print(f"Cancer Type: {cancer_type}\n Germline Event: {germline_event} \n Somatic Gene: {somatic_gene}")
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

    if germline_event not in df.columns or somatic_gene not in df.columns:       
        print(f"Combination {germline_event} - {somatic_gene} not found in DataFrame")
        return
    if df[germline_event].sum() == 0 or df[somatic_gene].sum() == 0:
        return

    df= df.dropna(subset=[germline_event,somatic_gene] + covariates)

    try:
        # Prepare the predictor (X) and response (y) variables
        X = df[[germline_event] + covariates]
        y = df[somatic_gene]
        
        # Normalize somatic_gene values: treat values > 1 as 1
        y = y.apply(lambda x: 1 if x > 1 else x)

        # Fit Logistic Regression Model
        model = FirthLogisticRegression()
        model.fit(X,y)
        
        # Extract p-value and odds ratio for the germline_gene predictor
        p_value = model.pvals_[0]
        odds_ratio = np.exp(model.coef_[0])
        
        # Calculate the 95% confidence interval for the odds ratio
        coef = model.coef_[0]
        std_err = model.bse_[0]

        # The confidence interval for the coefficient
        conf_int_coef = [coef - 1.96 * std_err, coef + 1.96 * std_err]

        # Convert the confidence interval from log odds to odds ratio
        conf_int_or = np.exp(conf_int_coef)

        # Append results to the results array, including confidence intervals
        return odds_ratio, p_value, conf_int_or[0], conf_int_or[1]

    except Exception as e:
        print(f"Error processing combination {germline_event} - {somatic_gene}: {e}")


def find_filtered_germline_event_frequency(df, cancer_type, germline_event, somatic_gene,germline_context,covariates=['male','pca_1','pca_2','pca_3','pca_4']):
    """
    Calculate the allele frequency for a specified column in a DataFrame.
    
    Args:
    df (pd.DataFrame): The input DataFrame.
    column_name (str): The column name for which to calculate allele frequency.
    
    Returns:
    float: The calculated allele frequency.
    """
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

    if germline_event not in df.columns or somatic_gene not in df.columns:       
        print(f"Combination {germline_event} - {somatic_gene} not found in DataFrame")
        return

    df= df.dropna(subset=[germline_event,somatic_gene] + covariates)
    
    # Get the column values, excluding NaNs
    column_values = df[germline_event].dropna()
    
    # Calculate the allele frequency
    total = None
    if germline_context == "coding":
        total = len(column_values)
    elif germline_context == "noncoding":
        total = (len(column_values) * 2)

    allele_frequency = column_values.sum() / total
    
    return allele_frequency, column_values.sum(), total

# Get a allele frequency for the cancer_type at large
def find_allele_frequency(df, cancer_type, germline_event):
    # Filter to Cancer type of interest
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

    # Make sure the snp is present in the df
    if germline_event not in df.columns:       
        print(f"Combination {germline_event} not found in DataFrame")
        return

    df= df.dropna(subset=[germline_event])
    
    # Get the column values, excluding NaNs
    column_values = df[germline_event].dropna()
    
    # Calculate the allele frequency
    allele_frequency = column_values.sum() / (len(column_values) * 2)
    
    return allele_frequency,column_values.sum(),(len(column_values) * 2)

# Get a mutation frequency for specific logistic regression (verify that all variables are there first)
def find_filtered_mutation_frequency(df, cancer_type, germline_event, somatic_gene, covariates=['male','pca_1','pca_2','pca_3','pca_4']):
    # Filter to Cancer Type of Interest
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

    # Verify that we have the data we want
    if germline_event not in df.columns or somatic_gene not in df.columns:       
        print(f"Combination {germline_event} - {somatic_gene} not found in DataFrame")
        return

    df= df.dropna(subset=[germline_event,somatic_gene] + covariates)
    
    # Get the column values, excluding NaNs
    column_values = df[somatic_gene].dropna()
    
    # Transform the column values: non-zero values become 1, zero values stay 0
    column_values = column_values.apply(lambda x: 1 if x != 0 else 0)
    
    # Calculate the allele frequency
    mutation_frequency = column_values.sum() / len(column_values)
    
    return mutation_frequency,column_values.sum(),len(column_values)

# Get a allele frequency for the cancer_type at large
def find_mutation_frequency(df, cancer_type, somatic_gene):
    # Filter to Cancer Type of Interest
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

    # Verify that we have the data we want
    if somatic_gene not in df.columns:       
        print(f"Combination {somatic_gene} not found in DataFrame")
        return

    df= df.dropna(subset=[somatic_gene])
    
    # Get the column values, excluding NaNs
    column_values = df[somatic_gene].dropna()
    
    # Transform the column values: non-zero values become 1, zero values stay 0
    column_values = column_values.apply(lambda x: 1 if x != 0 else 0)
    
    # Calculate the allele frequency
    mutation_frequency = column_values.sum() / len(column_values)
    
    return mutation_frequency,column_values.sum(),len(column_values)

def find_germline_event_frequency(df, cancer_type, germline_event, germline_context):
    """
    Calculate allele frequency for a germline event in a specific cancer type and context.

    Args:
        df (pd.DataFrame): Input DataFrame with cancer data.
        cancer_type (str): Cancer type of interest ('Pancancer' for all cancer types).
        germline_event (str): Germline event column to analyze.
        germline_context (str): Context of interest ('coding' or 'noncoding').

    Returns:
        tuple: (frequency, total germline events, total alleles considered)
    """
    # Filter to Cancer Type of Interest
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

    if df.empty:
        print("Error: find_germline_event_frequency - 0")
        return

    # Verify germline_event exists in the DataFrame
    if germline_event not in df.columns:
        print("Error: find_germline_event_frequency - 1")
        return

    # Remove rows with NaN in the germline_event column
    column_values = df[germline_event].dropna()

    if column_values.empty:
        print("Error: find_germline_event_frequency - 2")
        return

    # Calculate the allele frequency
    if germline_context == "coding":
        total = len(column_values)
    elif germline_context == "noncoding":
        total = len(column_values) * 2
    else:
        print("Error: find_germline_event_frequency - 3")
        return

    if total == 0:
        print("Error: find_germline_event_frequency - 4")
        return

    frequency = column_values.sum() / total

    return frequency, column_values.sum(), total


def merge_and_analyze_noncoding_coding(tsv1_path, tsv2_path, output_path='merged_with_combined_OR_and_p_values.tsv'):
    """
    This function merges two TSV files, performs inverse variance analysis on odds ratios (ORs),
    and combines p-values using Stouffer's Z-score method.

    Args:
    tsv1_path (str): File path to the HMF TSV file.
    tsv2_path (str): File path to the PROFILE TSV file.
    output_path (str): Path to save the output TSV file (default: 'merged_with_combined_OR_and_p_values.tsv').

    Returns:
    merged_df (DataFrame): The resulting DataFrame with combined OR and p-values.
    """

    # Step 1: Read the data into pandas DataFrames
    df1 = pd.read_csv(tsv1_path, sep='\t')
    df2 = pd.read_csv(tsv2_path, sep='\t')

    # Step 1: Replace 'Renal' with 'Kidney' in the 'cancer_type' column in both DataFrames
    df1['cancer_type'] = df1['cancer_type'].replace('Renal', 'Kidney')
    df2['cancer_type'] = df2['cancer_type'].replace('Renal', 'Kidney')

    # Step 2: Merge the DataFrames on the relevant columns, including 'criteria'
    merged_df = pd.merge(df1, df2, how='outer', 
                         on=['criteria', 'cancer_type', 'germline_risk_allele', 'germline_gene', 'germline_context', 'somatic_gene', 'somatic_context','relevant_cancer'],
                         suffixes=('_HMF', '_PROFILE'))

    # Step 3: Calculate variance for ORs based on confidence intervals (CI)
    merged_df['variance_HMF'] = ((np.log(merged_df['ci_OR_high_HMF']) - np.log(merged_df['ci_OR_low_HMF'])) / (2 * 1.96)) ** 2
    merged_df['variance_PROFILE'] = ((np.log(merged_df['ci_OR_high_PROFILE']) - np.log(merged_df['ci_OR_low_PROFILE'])) / (2 * 1.96)) ** 2

    # Step 4: Perform inverse variance weighting for OR
    merged_df['log_OR_HMF'] = np.log(merged_df['OR_HMF'])
    merged_df['log_OR_PROFILE'] = np.log(merged_df['OR_PROFILE'])

    # Calculate combined log OR using inverse variance weighting
    merged_df['log_OR_combined'] = (
        (merged_df['log_OR_HMF'] / merged_df['variance_HMF'] + merged_df['log_OR_PROFILE'] / merged_df['variance_PROFILE']) /
        (1 / merged_df['variance_HMF'] + 1 / merged_df['variance_PROFILE'])
    )

    # Convert combined log OR back to OR
    merged_df['OR_final'] = np.exp(merged_df['log_OR_combined'])

    # Step 5: Calculate the variance of the combined OR using inverse variance
    merged_df['variance_OR_combined'] = 1 / (1 / merged_df['variance_HMF'] + 1 / merged_df['variance_PROFILE'])

    # Step 6: Calculate Z-scores as OR_combined / variance_OR_combined
    merged_df['z_combined'] = np.log(merged_df['OR_final']) / np.sqrt(merged_df['variance_OR_combined'])

    # Step 7: Calculate p-value from Z-score
    merged_df['p_val_final'] = 2 * stats.norm.sf(np.abs(merged_df['z_combined']))

    # Step 8: Save the merged DataFrame with combined OR and p-values
    merged_df.to_csv(output_path, sep='\t', index=False)

    # Return the merged DataFrame
    return merged_df

def filter_by_pvalues(input_tsv, output_tsv = "filtered.tsv"):
    """
    This function filters the rows where (p_val_HMF or p_val_PROFILE < 0.05) 
    and p_val_combined < 0.05, and writes the filtered DataFrame to a new file.

    Args:
    input_tsv (str): Path to the input TSV file.
    output_tsv (str): Path to save the filtered output TSV file.
    
    Returns:
    filtered_df (DataFrame): The filtered DataFrame.
    """

    # Read the merged data
    df = pd.read_csv(input_tsv, sep='\t')

    # Filter based on the conditions
    filtered_df = df[
        ((df['p_val_HMF'] < 0.05) | (df['p_val_PROFILE'] < 0.05)) & 
        (df['p_val_combined'] < 0.05)
    ]

    # Save the filtered DataFrame to a new TSV file
    filtered_df.to_csv(output_tsv, sep='\t', index=False)

    # Return the filtered DataFrame
    return filtered_df

def analyze_data(convergence_table_path,genotype_table_path,germline_context,somatic_context,cohort,output_file):
  # Step 1: Read in Dataframe containing convergences and Dataframe containing patient level data  
  convergences_df = pd.read_csv(convergence_table_path,sep='\t')
  convergences_df = convergences_df[(convergences_df['germline_context'] == germline_context) & 
                                  (convergences_df['somatic_context'] == somatic_context)]
  patient_df = pd.read_csv(genotype_table_path, sep='\t')
  
  # Initialize an empty list to store results
  results = []

  # Step 2: Initialize covariates
  covariates = None
  if cohort == "HMF":
    covariates=['male','pca_1','pca_2','pca_3','age','TMB','tumorPurity','primary']
  elif cohort == "PROFILE":
    covariates=['male','pca_1','pca_2','pca_3','pca_4','age','TMB','tumorPurity','primary','late_stage']

  # Alternate slate of covariates for hormone sensitive cancers and pancancer
  covariates_hsc = [cov for cov in covariates if cov != 'male'] 
  covariates_pancancer = covariates + ['Breast_Diagnosis','Colorectal_Diagnosis','Kidney_Diagnosis','Lung_Diagnosis','Prostate_Diagnosis']
  
  # Step 3: Prepare Suffixes
  germline_suffix = None
  if germline_context == "coding":
    germline_suffix = "-g"
  else:
    germline_suffix = ""

  somatic_suffix = None
  if somatic_context == "coding":
    somatic_suffix = "-s"
  else:
    somatic_suffix = "-ncs"

  # Step 4: Loop through each row of the DataFrame and apply logistic regression with firth fallback
  for index, row in convergences_df.iterrows():
    germline_event = None

    # Extract necessary columns
    if germline_context == "coding":
      gwas_cancer_type = row[0]
      germline_gene = row[1]
      germline_context = row[2]
      somatic_gene = row[3]
      somatic_context = row[4]
      criteria = row[5]
      germline_event = germline_gene
    else:
      gwas_cancer_type = row[0]
      germline_snp = row[1]
      germline_gene = row[2]
      germline_context = row[3]
      somatic_gene = row[4]
      somatic_context = row[5]
      criteria = row[6]
      germline_event = germline_snp
    
    for cancer_type in ["Breast","Colorectal","Prostate","Lung","Kidney","Pancancer"]:
      relevant_cancer = 0
      if gwas_cancer_type.lower() == cancer_type.lower():
        relevant_cancer = 1
      elif cancer_type == "Pancancer":
        relevant_cancer = 2

      # Apply the functions to get metrics
      if cancer_type == "Breast" or cancer_type == "Prostate":
        regression_output = logistic_regression_with_fallback(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix,covariates=covariates_hsc)
        filtered_germline_output = find_filtered_germline_event_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix,germline_context, covariates=covariates_hsc)
        filtered_somatic_output = find_filtered_mutation_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, covariates=covariates_hsc)
      elif cancer_type == "Pancancer":
        regression_output = logistic_regression_with_fallback(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix,covariates=covariates_pancancer)
        filtered_germline_output = find_filtered_germline_event_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, germline_context, covariates=covariates_pancancer)
        filtered_somatic_output = find_filtered_mutation_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, covariates=covariates_pancancer)       
      else:
        regression_output = logistic_regression_with_fallback(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix,covariates=covariates)
        filtered_germline_output = find_filtered_germline_event_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, germline_context, covariates=covariates)
        filtered_somatic_output = find_filtered_mutation_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, covariates=covariates)
              
      germline_output = find_germline_event_frequency(patient_df, cancer_type, germline_event + germline_suffix,germline_context)
      somatic_output = find_mutation_frequency(patient_df, cancer_type, somatic_gene + somatic_suffix)
          
      # Assign values using ternary operators for cleaner code
      OR, p_val, ci_OR_low, ci_OR_high = regression_output if regression_output else (pd.NA, pd.NA, pd.NA, pd.NA)

      filtered_germline_plp_frequency, filtered_germline_plp_count, filtered_germline_sample_size = (
        filtered_germline_output if filtered_germline_output else (pd.NA, pd.NA, pd.NA)
      )

      germline_plp_frequency, germline_plp_count, germline_sample_size = (
        germline_output if germline_output else (pd.NA, pd.NA, pd.NA)
      )

      somatic_plp_frequency, somatic_plp_count, somatic_sample_size = (
        somatic_output if somatic_output else (pd.NA, pd.NA, pd.NA)
      )

      filtered_somatic_plp_frequency, filtered_somatic_plp_count, filtered_somatic_sample_size = (
        filtered_somatic_output if filtered_somatic_output else (pd.NA, pd.NA, pd.NA)
      )

      if germline_context == "coding":
        # Append the results to the list
        results.append([
            cancer_type, germline_gene, germline_context, 
            somatic_gene, somatic_context, criteria,
            germline_plp_frequency, germline_plp_count, germline_sample_size,
            filtered_germline_plp_frequency, filtered_germline_plp_count, filtered_germline_sample_size,
            somatic_plp_frequency, somatic_plp_count, somatic_sample_size,
            filtered_somatic_plp_frequency, filtered_somatic_plp_count, filtered_somatic_sample_size,
            OR, p_val, ci_OR_low, ci_OR_high, relevant_cancer
        ])
      else:
        # Append the results to the list
        results.append([
            cancer_type, germline_snp, germline_gene, germline_context, 
            somatic_gene, somatic_context, criteria,
            germline_plp_frequency, germline_plp_count, germline_sample_size,
            filtered_germline_plp_frequency, filtered_germline_plp_count, filtered_germline_sample_size,
            somatic_plp_frequency, somatic_plp_count, somatic_sample_size,
            filtered_somatic_plp_frequency, filtered_somatic_plp_count, filtered_somatic_sample_size,
            OR, p_val, ci_OR_low, ci_OR_high, relevant_cancer
        ])          
  
  if germline_context == "noncoding" and somatic_context == "coding":
    cohort_suffix = cohort
  else:
    cohort_suffix = "final"

  if germline_context == "coding":
    # Convert the results to a DataFrame
    output_df = pd.DataFrame(results, columns=[
        'cancer_type', 'germline_gene', 'germline_context', 
        'somatic_gene', 'somatic_context', 'criteria',
        f'germline_plp_frequency_{cohort}', f'germline_plp_count_{cohort}', f'germline_sample_size_{cohort}',
        f'filtered_germline_plp_frequency_{cohort}', f'filtered_germline_plp_count_{cohort}', f'filtered_germline_sample_size_{cohort}',
        f'somatic_plp_frequency_{cohort}', f'somatic_plp_count_{cohort}', f'somatic_sample_size_{cohort}',
        f'filtered_somatic_plp_frequency_{cohort}', f'filtered_somatic_plp_count_{cohort}', f'filtered_somatic_sample_size_{cohort}',
        f'OR_{cohort_suffix}', f'p_val_{cohort_suffix}', f'ci_OR_low_{cohort_suffix}', f'ci_OR_high_{cohort_suffix}', 'relevant_cancer'
    ])
  else:
    # Convert the results to a DataFrame
    output_df = pd.DataFrame(results, columns=[
        'cancer_type', 'germline_risk_allele','germline_gene', 'germline_context', 
        'somatic_gene', 'somatic_context', 'criteria',
        f'germline_plp_frequency_{cohort}', f'germline_plp_count_{cohort}', f'germline_sample_size_{cohort}',
        f'filtered_germline_plp_frequency_{cohort}', f'filtered_germline_plp_count_{cohort}', f'filtered_germline_sample_size_{cohort}',
        f'somatic_plp_frequency_{cohort}', f'somatic_plp_count_{cohort}', f'somatic_sample_size_{cohort}',
        f'filtered_somatic_plp_frequency_{cohort}', f'filtered_somatic_plp_count_{cohort}', f'filtered_somatic_sample_size_{cohort}',
        f'OR_{cohort_suffix}', f'p_val_{cohort_suffix}', f'ci_OR_low_{cohort_suffix}', f'ci_OR_high_{cohort_suffix}', 'relevant_cancer'
    ])

  # Eliminate duplicate rows based on all columns
  output_df = output_df.drop_duplicates()

  # Write the DataFrame to a TSV file
  output_df.to_csv(output_file, sep='\t', index=False)

def process_and_merge_profile_tables(variant_table_pathway, gene_table_pathway):
    """
    Process variant and gene tables by filtering columns based on suffixes, merging them,
    and adjusting '-s' columns based on '-g' column values.

    Parameters:
    - variant_table: DataFrame containing variant data.
    - gene_table: DataFrame containing gene data.

    Returns:
    - DataFrame: The processed and merged table.
    """

    variant_table = pd.read_csv(variant_table_pathway,sep='\t')
    gene_table = pd.read_csv(gene_table_pathway,sep='\t')

    # Step 1: Process column names
    # Remove '-s' and '-g' suffixes for comparison
    variant_columns = [col for col in variant_table.columns if col != 'patient_id' and (col.endswith('-s') or col.endswith('-g'))]
    gene_columns = [col for col in gene_table.columns if col != 'patient_id' and (col.endswith('-s') or col.endswith('-g'))]

    variant_base_columns = [col.rsplit('-', 1)[0] for col in variant_columns]
    gene_base_columns = [col.rsplit('-', 1)[0] for col in gene_columns]

    # Step 2: Find common base columns
    common_columns = list(set(variant_base_columns).intersection(set(gene_base_columns)))

    # Step 3: Create filtered_gene_columns and filter gene_table
    filtered_gene_columns = [col + '-g' for col in common_columns] + ['patient_id']
    filtered_gene_table = gene_table[filtered_gene_columns]

    # Step 4: Merge filtered_gene_table with variant_table
    merged_table = variant_table.merge(filtered_gene_table, on='patient_id', how='inner')

    # Step 5: Adjust '-s' columns based on '-g' columns
    conversion_stats = []  # To store stats for each column
    for col in common_columns:
        s_col = col + '-s'
        g_col = col + '-g'
        if s_col in merged_table.columns and g_col in merged_table.columns:
            # Get initial count of values > 0 in the s_col
            initial_s_col_count = (merged_table[s_col] > 0).sum()

            # Perform the conversion
            merged_table.loc[(merged_table[g_col] == 1) & (merged_table[s_col] > 0), s_col] = 0

            # Get the count of values converted to 0
            converted_count = initial_s_col_count - (merged_table[s_col] > 0).sum()

            # Calculate percentage of converted values
            if initial_s_col_count > 0:
                conversion_percentage = (converted_count / initial_s_col_count) * 100
            else:
                conversion_percentage = 0.0

            # Append stats to the list
            conversion_stats.append({
                "s_col": s_col,
                "g_col": g_col,
                "converted_count": converted_count,
                "initial_s_col_count": initial_s_col_count,
                "conversion_percentage": conversion_percentage
            })

    # Print stats for each column
    for stats in conversion_stats:
        print(f"Column: {stats['s_col']} | "
              f"Converted: {stats['converted_count']} out of {stats['initial_s_col_count']} "
              f"({stats['conversion_percentage']:.2f}%)")


    merged_table.to_csv("profile_table.tsv",sep='\t',index=False)
    return

