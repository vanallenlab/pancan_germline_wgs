import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact
import numpy as np
#from firthlogist import FirthLogisticRegression

def perform_logistic_regression(df, combination, results_array):
    germline_gene, somatic_gene = combination
    
    # Check if both columns exist in the DataFrame
    if germline_gene not in df.columns:
        print(f"germline_gene {germline_gene} not found in DataFrame columns")
        return
    if somatic_gene not in df.columns:
        print(f"somatic_gene {somatic_gene} not found in DataFrame columns")
        return

    # Filter out rows with NA values in germline or somatic gene columns
    df_filtered = df.dropna(subset=[germline_gene, somatic_gene])
    
    if df_filtered[germline_gene].sum() == 0:
        #print(f"No germline_gene: {germline_gene}")
        return
    if df_filtered[somatic_gene].sum() == 0:
        #print(f"No somatic_gene: {somatic_gene}")
        return
    
    try:
        # Prepare the predictor (X) and response (y) variables
        X = df_filtered[[germline_gene]]
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
        p_value = result.pvalues[germline_gene]
        odds_ratio = np.exp(result.params[germline_gene])
        
        # Calculate the 95% confidence interval for the odds ratio
        coef = result.params[germline_gene]
        std_err = result.bse[germline_gene]

        # The confidence interval for the coefficient
        conf_int_coef = [coef - 1.96 * std_err, coef + 1.96 * std_err]

        # Convert the confidence interval from log odds to odds ratio
        conf_int_or = np.exp(conf_int_coef)

        # Append results to the results array, including confidence intervals
        results_array.append((combination, p_value, odds_ratio, conf_int_or[0], conf_int_or[1], af, ac, an, mf, mc, mn))

    except Exception as e:
        print(f"Error processing combination {combination}: {e}")

def get_unique_combinations_nc_c(tsv_file,germline_context,somatic_context):
    unique_combinations = set()
    if germline_context == "coding" and somatic_context == "coding":
        with open(tsv_file, 'r') as file:
            # Skip the header line
            next(file)
        
            for line in file:
                line = line.strip().split('\t')
                germline_gene = line[1]+'-g'
                germline_context = line[2]
                somatic_mutation = line[3]+'-s'
                somatic_context = line[4]
                criteria = line[5]
                combination = (germline_gene, somatic_mutation)
                unique_combinations.add((combination,criteria))
    if germline_context == "noncoding" and somatic_context == "coding":
        with open(tsv_file, 'r') as file:
            # Skip the header line
            next(file)
        
            for line in file:
                line = line.strip().split('\t')
                germline_risk_allele = line[1]
                somatic_mutation = line[4]+'-s'
                criteria = line[5]
                combination = (germline_risk_allele, somatic_mutation)
                unique_combinations.add((combination,criteria))
    elif germline_context == "noncoding" and somatic_context == "noncoding":
    	pass
    elif germline_context == "coding" and somatic_context == "noncoding":
    	pass
    else:
    	raise ValueError("An error occurred: Germline and Somatic Context have to be either 'coding' or 'noncoding'. ")

    return unique_combinations


def fishers_exact(df, cancer_type, germline_gene,somatic_gene):
    if cancer_type != "pancancer":
        df = df[df['cancer_type'] == cancer_type]
    germline_gene_full_name = germline_gene + '-g'
    somatic_gene_full_name = somatic_gene + '-s'

    if germline_gene_full_name not in df.columns or somatic_gene_full_name not in df.columns:       
        print(f"Combination {combination} not found in DataFrame")
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
        log_or = np.log(odds_ratio)
        table_array = np.array(table)
        
        # Using a standard error approximation for the log odds ratio
        a, b, c, d = table_array.ravel()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        
        # Compute 95% confidence interval for log(odds ratio)
        log_ci_low = log_or - 1.96 * se
        log_ci_high = log_or + 1.96 * se

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

def find_mutation_frequency(df, column_name):
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
    
    # Transform the column values: non-zero values become 1, zero values stay 0
    column_values = column_values.apply(lambda x: 1 if x != 0 else 0)
    
    # Calculate the allele frequency
    mutation_frequency = column_values.sum() / (len(column_values) * 2)
    
    return mutation_frequency,column_values.sum(),len(column_values)

