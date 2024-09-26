import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact
import numpy as np
from firthlogist import FirthLogisticRegression

def firth_logistic_regression(df, cancer_type, germline_event, somatic_gene):
    print(f"Cancer Type: {cancer_type}\n Germline Event: {germline_event} \n Somatic Gene: {somatic_gene}")
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]

    df= df.dropna(subset=['male', 'pca_1', 'pca_2', 'pca_3', 'pca_4', 'stage'])

    if germline_event not in df.columns or somatic_gene not in df.columns:       
        print(f"Combination {germline_event} - {somatic_gene} not found in DataFrame")
        return
    if df[germline_event].sum() == 0 or df[somatic_gene].sum() == 0:
        return

    # Check if both columns exist in the DataFrame
    if germline_event not in df.columns:
        print(f"germline_gene {germline_event} not found in DataFrame columns")
        return
    if somatic_gene not in df.columns:
        print(f"somatic_gene {somatic_gene} not found in DataFrame columns")
        return
    
    try:
        # Prepare the predictor (X) and response (y) variables
        X = df[[germline_event,'male','pca_1','pca_2','pca_3','pca_4','stage']]
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


def find_allele_frequency(df, column_name):
    """
    Calculate the allele frequency for a specified column in a DataFrame.
    
    Args:
    df (pd.DataFrame): The input DataFrame.
    column_name (str): The column name for which to calculate allele frequency.
    
    Returns:
    float: The calculated allele frequency.
    """
    # Check if the column exists in the DataFrame
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in DataFrame.")
    
    # Get the column values, excluding NaNs
    column_values = df[column_name].dropna()
    
    # Calculate the allele frequency
    allele_frequency = column_values.sum() / (len(column_values) * 2)
    
    return allele_frequency,column_values.sum(),(len(column_values) * 2)

def find_mutation_frequency(df, cancer_type, column_name):
    """
    Calculate the allele frequency for a specified column in a DataFrame.
    
    Args:
    df (pd.DataFrame): The input DataFrame.
    column_name (str): The column name for which to calculate allele frequency.
    
    Returns:
    float: The calculated allele frequency.
    """
    # Check if the column exists in the DataFrame
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]
    df= df.dropna(subset=['male', 'pca_1', 'pca_2', 'pca_3', 'pca_4', 'stage'])

    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in DataFrame.")
    
    # Get the column values, excluding NaNs
    column_values = df[column_name].dropna()
    
    # Transform the column values: non-zero values become 1, zero values stay 0
    column_values = column_values.apply(lambda x: 1 if x != 0 else 0)
    
    # Calculate the allele frequency
    mutation_frequency = column_values.sum() / len(column_values)
    
    return mutation_frequency,column_values.sum(),len(column_values)

