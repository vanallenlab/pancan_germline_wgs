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
                         on=['criteria', 'cancer_type', 'germline_risk_allele', 'germline_gene', 'germline_context', 'somatic_gene', 'somatic_context'],
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
    merged_df['OR_combined'] = np.exp(merged_df['log_OR_combined'])

    # Step 5: Calculate the variance of the combined OR using inverse variance
    merged_df['variance_OR_combined'] = 1 / (1 / merged_df['variance_HMF'] + 1 / merged_df['variance_PROFILE'])

    # Step 6: Calculate Z-scores as OR_combined / variance_OR_combined
    merged_df['z_combined'] = np.log(merged_df['OR_combined']) / np.sqrt(merged_df['variance_OR_combined'])

    # Step 7: Calculate p-value from Z-score
    merged_df['p_val_combined'] = 2 * stats.norm.sf(np.abs(merged_df['z_combined']))

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

