import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact
import numpy as np
from firthlogist import FirthLogisticRegression
import scipy.stats as stats

import statsmodels.api as sm
import numpy as np

# Assuming FirthLogisticRegression is already imported or defined elsewhere
# You may need to implement this if it's not available.

def logistic_regression_with_fallback(df, cancer_type, germline_event, somatic_gene, covariates=['male', 'pca_1', 'pca_2', 'pca_3', 'pca_4']):
    print(f"Cancer Type: {cancer_type}\n Germline Event: {germline_event} \n Somatic Gene: {somatic_gene}")
    
    # Filter based on cancer_type if it's not Pancancer
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

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


def find_filtered_allele_frequency(df, cancer_type, germline_event, somatic_gene,covariates=['male','pca_1','pca_2','pca_3','pca_4']):
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
    allele_frequency = column_values.sum() / (len(column_values) * 2)
    
    return allele_frequency,column_values.sum(),(len(column_values) * 2)

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
  convergences_df = convergences_df[convergences_df['germline_context'] == germline_context and convergences_df['somatic_context'] == somatic_context]
  genotype_df = pd.read_csv(genotype_table_path, sep='\t')
  
  # Initialize an empty list to store results
  results = []

  # Step 2: Initialize covariates
  covariates = None
  if cohort == "HMF":
    covariates=['male','pca_1','pca_2','pca_3','pca_4']
  elif cohort == "PROFILE":
    covariates=['male','pca_1','pca_2','pca_3','pca_4','late_stage']

  # Alternate slate of covariates for hormone sensitive cancers
  covariates_hsc = [cov for cov in covariates if cov != 'male'] 
  
  # Step 3: Prepare Suffixes
  germline_suffix = None
  if germline_context == "coding":
    germline_suffix = "-g"

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
        regression_output = gsc_util.logistic_regression_with_fallback(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix,covariates=covariates_hsc)
        filtered_germline_output = gsc_util.find_filtered_allele_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, covariates=covariates_hsc)
        filtered_somatic_output = gsc_util.find_filtered_mutation_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, covariates=covariates_hsc)
      else:
        regression_output = gsc_util.logistic_regression_with_fallback(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix,covariates=covariates)
        filtered_germline_output = gsc_util.find_filtered_allele_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, covariates=covariates)
        filtered_somatic_output = gsc_util.find_filtered_mutation_frequency(patient_df, cancer_type, germline_event + germline_suffix, somatic_gene + somatic_suffix, covariates=covariates)
              
        germline_output = gsc_util.find_mutation_frequency(patient_df, cancer_type, germline_event + germline_suffix)
        somatic_output = gsc_util.find_mutation_frequency(patient_df, cancer_type, somatic_gene + somatic_suffix)
          
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
          cancer_type, germline_snp, germline_gene, germline_context, 
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
          cancer_type, germline_gene, germline_context, 
          somatic_gene, somatic_context, criteria,
          germline_plp_frequency, germline_plp_count, germline_sample_size,
          filtered_germline_plp_frequency, filtered_germline_plp_count, filtered_germline_sample_size,
          somatic_plp_frequency, somatic_plp_count, somatic_sample_size,
          filtered_somatic_plp_frequency, filtered_somatic_plp_count, filtered_somatic_sample_size,
          OR, p_val, ci_OR_low, ci_OR_high, relevant_cancer
        ])          

  if germline_context == "coding":
    # Convert the results to a DataFrame
    output_df = pd.DataFrame(results, columns=[
        'cancer_type', 'germline_risk_allele','germline_gene', 'germline_context', 
        'somatic_gene', 'somatic_context', 'criteria',
        f'germline_plp_frequency_{cohort}', f'germline_plp_count_{cohort}', f'germline_sample_size_{cohort}',
        f'filtered_germline_plp_frequency_{cohort}', f'filtered_germline_plp_count_{cohort}', f'filtered_germline_sample_size_{cohort}',
        f'somatic_plp_frequency_{cohort}', f'somatic_plp_count_{cohort}', f'somatic_sample_size_{cohort}',
        f'filtered_somatic_plp_frequency_{cohort}', f'filtered_somatic_plp_count_{cohort}', f'filtered_somatic_sample_size_{cohort}',
        'OR_final', 'p_val_final', 'ci_OR_low_final', 'ci_OR_high_final', 'relevant_cancer'
    ])
  else:
    # Convert the results to a DataFrame
    output_df = pd.DataFrame(results, columns=[
        'cancer_type', 'germline_gene', 'germline_context', 
        'somatic_gene', 'somatic_context', 'criteria',
        f'germline_plp_frequency_{cohort}', f'germline_plp_count_{cohort}', f'germline_sample_size_{cohort}',
        f'filtered_germline_plp_frequency_{cohort}', f'filtered_germline_plp_count_{cohort}', f'filtered_germline_sample_size_{cohort}',
        f'somatic_plp_frequency_{cohort}', f'somatic_plp_count_{cohort}', f'somatic_sample_size_{cohort}',
        f'filtered_somatic_plp_frequency_{cohort}', f'filtered_somatic_plp_count_{cohort}', f'filtered_somatic_sample_size_{cohort}',
        'OR_final', 'p_val_final', 'ci_OR_low_final', 'ci_OR_high_final', 'relevant_cancer'
    ])

  # Eliminate duplicate rows based on all columns
  output_df = output_df.drop_duplicates()

  # Write the DataFrame to a TSV file
  output_df.to_csv(output_file, sep='\t', index=False)
