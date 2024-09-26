import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact
import numpy as np
from firthlogist import FirthLogisticRegression

def firth_logistic_regression(df, cancer_type, germline_event, somatic_gene):
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]
    # Check if both columns exist in the DataFrame
    if germline_event not in df.columns:
        print(f"germline_gene {germline_event} not found in DataFrame columns")
        return
    if somatic_gene not in df.columns:
        print(f"somatic_gene {somatic_gene} not found in DataFrame columns")
        return

    # Filter out rows with NA values in germline or somatic gene columns
    df_filtered = df.dropna(subset=[germline_event, somatic_gene])
    
    if df_filtered[germline_event].sum() == 0 or df_filtered[somatic_gene].sum() == 0:
        #print(f"No germline_gene: {germline_gene}")
        return
    
    try:
        # Prepare the predictor (X) and response (y) variables
        X = df_filtered[[germline_event,male,pca_1,pca_2,pca_3,pca_4,stage]]
        y = df_filtered[somatic_gene]
        
        # Normalize somatic_gene values: treat values > 1 as 1
        y = y.apply(lambda x: 1 if x > 1 else x)

        # Add a constant to the predictor variables
        X = sm.add_constant(X)

        # Fit the logistic regression model with disp=False to suppress output
        #model = sm.Logit(y, X)
        #result = model.fit(disp=False)
        model = FirthLogisticRegression()
        model.fit(X,y)
        
        # Extract p-value and odds ratio for the germline_gene predictor
        p_value = result.pvalues[germline_event]
        odds_ratio = np.exp(result.params[germline_event])
        
        # Calculate the 95% confidence interval for the odds ratio
        coef = result.params[germline_event]
        std_err = result.bse[germline_event]

        # The confidence interval for the coefficient
        conf_int_coef = [coef - 1.96 * std_err, coef + 1.96 * std_err]

        # Convert the confidence interval from log odds to odds ratio
        conf_int_or = np.exp(conf_int_coef)

        # Append results to the results array, including confidence intervals
        return p_value, odds_ratio, conf_int_or[0], conf_int_or[1]

    except Exception as e:
        print(f"Error processing combination {germline_event} - {somatic_gene}: {e}")


def fishers_exact(df, cancer_type, germline_gene,somatic_gene):
    if cancer_type != "Pancancer":
        df = df[df['cancer_type'] == cancer_type]
    germline_gene_full_name = germline_gene + '-g'
    somatic_gene_full_name = somatic_gene + '-s'

    if germline_gene_full_name not in df.columns or somatic_gene_full_name not in df.columns:       
        print(f"Combination {germline_gene} - {somatic_gene} not found in DataFrame")
        return
    if df[germline_gene_full_name].sum() == 0 or df[somatic_gene_full_name].sum() == 0:
        return

    try:
        # Create a 2x2 contingency table manually
        table = [
            [
                df[(df[germline_gene_full_name] == 0) & (df[somatic_gene_full_name] == 0)].shape[0],
                df[(df[germline_gene_full_name] == 0) & (df[somatic_gene_full_name] == 1)].shape[0]
            ],
            [
                df[(df[germline_gene_full_name] == 1) & (df[somatic_gene_full_name] == 0)].shape[0],
                df[(df[germline_gene_full_name] == 1) & (df[somatic_gene_full_name] == 1)].shape[0]
            ]
        ]

        
        # Perform Fisher's exact test
        odds_ratio, p_value = fisher_exact(table)

        # Check if odds_ratio is zero
        if odds_ratio == 0:
            log_or = -np.inf  # Use negative infinity for log odds ratio
        else:
            log_or = np.log(odds_ratio)

        # Convert the table to a NumPy array
        table_array = np.array(table)

        # Standard error approximation for the log odds ratio
        a, b, c, d = table_array.ravel()
        if a == 0 or b == 0 or c == 0 or d == 0:
            epsilon = 0.01
            a += epsilon
            b += epsilon
            c += epsilon
            d += epsilon

        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

        # Compute 95% confidence interval for log(odds ratio)
        log_ci_low = log_or - 1.96 * se
        log_ci_high = log_or + 1.96 * se

        # If log_or is -inf, ci_low will be very negative and should be handled
        if log_or == -np.inf:
            ci_low = 0  # You might set the lower bound of the CI to 0
            ci_high = 0  # Or set it to something reasonable, e.g., 0
        else:
            ci_low = np.exp(log_ci_low)
            ci_high = np.exp(log_ci_high)

        return odds_ratio, p_value, ci_low, ci_high
        
    except Exception as e:
        print(f"Error occurred for combination {germline_gene} & {somatic_gene}: {str(e)}")
        return

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
        
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in DataFrame.")
    
    # Get the column values, excluding NaNs
    column_values = df[column_name].dropna()
    
    # Transform the column values: non-zero values become 1, zero values stay 0
    column_values = column_values.apply(lambda x: 1 if x != 0 else 0)
    
    # Calculate the allele frequency
    mutation_frequency = column_values.sum() / len(column_values)
    
    return mutation_frequency,column_values.sum(),len(column_values)

