import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from statsmodels.stats.multitest import fdrcorrection

def plot_germline_frequencies(data, 
                              hmf_col='germline_plp_frequency_HMF', 
                              profile_col='germline_plp_frequency_PROFILE', 
                              save_path="germline_snp_frequencies.png"):

    # Define the cancer types of interest
    cancer_types = ['Breast', 'Prostate', 'Colorectal', 'Lung', 'Kidney', 'Pancancer']
    
    # Drop rows where either the HMF or PROFILE frequency is missing (NaN)
    data_clean = data.dropna(subset=[hmf_col, profile_col])
    data = data_clean

    # Filtering to keep only the specified columns
    data = data[['germline_risk_allele', 'cancer_type', hmf_col, profile_col]]

    # Removing duplicates
    data = data.drop_duplicates()

    # Create a figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()  # Flatten the array for easy indexing

    # Iterate through each cancer type and create a subplot
    for idx, cancer_type in enumerate(cancer_types):

        # Filter the DataFrame for the current cancer type
        filtered_data = data[data['cancer_type'] == cancer_type]
        print(cancer_type)
        # Extract HMF and PROFILE frequencies
        x = filtered_data[hmf_col].values.reshape(-1, 1)
        y = filtered_data[profile_col].values

        # Fit a linear regression model
        model = LinearRegression()
        model.fit(x, y)
        y_pred = model.predict(x)
        r_squared = r2_score(y, y_pred)

        # Plot the data
        axes[idx].scatter(x, y, alpha=0.5)
        axes[idx].plot(x, y_pred, color='blue', label='Best Fit', linestyle='--')
        
        # Add y=x line
        axes[idx].plot(x, x, color='black', linestyle='--', label='y = x')

        # Set titles and labels
        axes[idx].set_title(f'{cancer_type} (R² = {r_squared:.2f})')
        axes[idx].set_xlabel('Germline Allele Frequency HMF')
        axes[idx].set_ylabel('Germline Allele Frequency PROFILE')
        axes[idx].legend()
        axes[idx].set_xlim(0, max(x.max(), y.max()))
        axes[idx].set_ylim(0, max(x.max(), y.max()))

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_somatic_mutation_frequencies(df, 
                                        x_col="somatic_plp_frequency_HMF", 
                                        y_col="somatic_plp_frequency_PROFILE", 
                                        save_path="somatic_mutation_frequencies.png"):
    # Specify the cancer types to include
    cancer_types = ["Breast", "Prostate", "Colorectal", "Lung", "Kidney", "Pancancer"]

    # Drop rows where either the HMF or PROFILE frequency is missing (NaN)
    df_clean = df.dropna(subset=[x_col, y_col])
    df = df_clean

    # Filtering to keep only the specified columns
    df = df[['somatic_gene', 'cancer_type', x_col, y_col]]

    # Removing duplicates
    df = df.drop_duplicates()

    # Set up the figure and axes
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for i, cancer_type in enumerate(cancer_types):
        # Filter data for the current cancer type
        subset = df[df['cancer_type'] == cancer_type]
        
        # Prepare data for regression
        X = subset[x_col].values.reshape(-1, 1)
        y = subset[y_col].values

        # Fit linear regression
        model = LinearRegression()
        model.fit(X, y)
        
        # Get predictions and R² value
        y_pred = model.predict(X)
        r2 = r2_score(y, y_pred)
        
        # Plot the data
        axes[i].scatter(X, y, alpha=0.6, label=cancer_type)
        axes[i].plot(X, y_pred, linestyle='--', color='red', label='Best Fit')
        axes[i].plot([X.min(), X.max()], [X.min(), X.max()], 'k--', label='y=x')  # dashed line for y=x
        axes[i].set_title(f"{cancer_type}: R² = {r2:.2f}")
        axes[i].set_xlabel(x_col)
        axes[i].set_ylabel(y_col)
        axes[i].legend()
        
    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_volcano(df, 
                 p_col="p_val_combined", 
                 or_col="OR_combined", 
                 criteria_col="criteria", 
                 cancer_type_col="cancer_type",
                 germline_context="germline_risk_allele",
                 figure_title="Insert Title Here",
                 save_path="insert_path_here.png"):

    # Drop rows where either the HMF or PROFILE frequency is missing (NaN)
    df_clean = df.dropna(subset=[p_col, or_col])
    df = df_clean
    
    # Calculate -log10(p-value) and log2(OR)
    df['log10_p'] = -np.log10(df[p_col])
    df['log2_OR'] = np.log2(df[or_col])
    
    # Set up point shapes based on 'criteria' and colors based on 'cancer_type'
    shapes = {'known_ppi': 'o', 'same_gene': 's', 'protein_complex':'h','ligand_receptor':'^',
            'known_ppi;ligand_receptor':'*','known_ppi;protein_complex':'d'}  # Example shapes for criteria
    colors = {'Breast': 'red', 'Prostate': 'blue', 'Colorectal': 'green', 
              'Lung': 'purple', 'Kidney': 'orange', 'Pancancer': 'black'}
    
    # Create figure and axes
    plt.figure(figsize=(10, 8))
    
    # Plot points with different shapes and colors
    for criteria_val, shape in shapes.items():
        for cancer_type, color in colors.items():
            subset = df[(df[criteria_col] == criteria_val) & (df[cancer_type_col] == cancer_type)]
            plt.scatter(subset['log2_OR'], subset['log10_p'], 
                        marker=shape, color=color, label=f"{criteria_val}-{cancer_type}", alpha=0.7)
    
    # Add Bonferroni and nominal significance lines
    bonferroni_threshold = -np.log10(0.05 / len(df))  # Bonferroni threshold
    nominal_threshold = -np.log10(0.05)  # Nominal significance threshold
    rejected, pvals_corrected = fdrcorrection(df[p_col], alpha=0.05, method='indep', is_sorted=False)

    # We calculate the highest corrected p-value that is still below the alpha threshold (0.05)
    # Check if any values are rejected and handle the case where none are
    if rejected.any():  # If any p-values are below the FDR threshold
        fdr_threshold = -np.log10(pvals_corrected[rejected].max())
        plt.axhline(y=fdr_threshold, color='green', linestyle='--', label=f'FDR Threshold: {-np.log10(fdr_threshold):.2f}')
    else:
        fdr_threshold = -np.log10(0.05)  # Or set it to a default value, e.g., np.nan

    # Optionally, print a warning if no p-values passed the FDR threshold
    if fdr_threshold is None:
        print("Warning: No p-values passed the FDR threshold.")
        
    plt.axhline(y=bonferroni_threshold, color='red', linestyle='--', label=f"Bonferroni Significance (p = {round(0.05 / len(df),7)})")
    plt.axhline(y=nominal_threshold, color='blue', linestyle='--', label='Nominal Significance (p=0.05)')
    
    # Label points above Bonferroni threshold
    significant_points = df[df['log10_p'] > nominal_threshold]
    for _, row in significant_points.iterrows():
        if pd.notna(row[germline_context]):  # Check if germline_context is NOT NaN
            label_text = f"({row[germline_context]} ,{row['somatic_gene']})"
        else:
            label_text = f"({row['germline_gene']} ,{row['somatic_gene']})"
        plt.text(row['log2_OR'], row['log10_p'], label_text, fontsize=8, ha='right')

    # Set plot labels and title
    plt.xlabel('log2(OR)')
    plt.ylabel('-log10(p)')
    plt.title(f'{figure_title}')
    
    # Display legend
    # Add a custom legend for cancer_type (colors)
    color_handles = [plt.Line2D([0], [0], marker='o', color='w', label=cancer, 
                                markerfacecolor=color, markersize=10) 
                     for cancer, color in colors.items()]
    plt.legend(handles=color_handles, title="Cancer Type", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add a separate custom legend for criteria (shapes) below the cancer_type legend
    shape_handles = [plt.Line2D([0], [0], marker=shape, color='k', label=criteria, markersize=10)
                     for criteria, shape in shapes.items()]
    plt.legend(handles=shape_handles, title="Criteria", bbox_to_anchor=(1.05, 0.75), loc='upper left')

    # Set x-axis limit and add vertical dashed line at x=0
    plt.xlim(left=-20)  # Limit the left side of the x-axis to -20
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.3)  # Light dashed line

    # Save plot to file
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def plot_p_values(data, 
                  hmf_col='p_val_HMF', 
                  profile_col='p_val_PROFILE', 
                  save_path="p_val_comparison.png"):
    """
    Function to plot -log10 transformed p-values for HMF and PROFILE across different cancer types.
    """
    # Define the cancer types of interest
    cancer_types = ['Breast', 'Prostate', 'Colorectal', 'Lung', 'Kidney', 'Pancancer']

    # Drop rows where either the HMF or PROFILE p-value is missing (NaN)
    data_clean = data.dropna(subset=[hmf_col, profile_col])
    data_clean = data_clean[(data_clean[hmf_col] > 0) & (data_clean[profile_col] > 0)]
    data = data_clean

    # Transform p-values to -log10(p)
    data['log_p_val_HMF'] = -np.log10(data[hmf_col])
    data['log_p_val_PROFILE'] = -np.log10(data[profile_col])

    # Drop rows with non-finite values after transformation
    data = data_clean[np.isfinite(data['log_p_val_HMF']) & np.isfinite(data['log_p_val_PROFILE'])]

    # Create a figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()  # Flatten the array for easy indexing

    # Iterate through each cancer type and create a subplot
    for idx, cancer_type in enumerate(cancer_types):
        
        # Filter the DataFrame for the current cancer type
        filtered_data = data[data['cancer_type'] == cancer_type]

        # Extract transformed HMF and PROFILE p-values
        x = filtered_data['log_p_val_HMF'].values.reshape(-1, 1)
        y = filtered_data['log_p_val_PROFILE'].values

        # Fit a linear regression model
        model = LinearRegression()
        model.fit(x, y)
        y_pred = model.predict(x)
        r_squared = r2_score(y, y_pred)

        # Plot the data
        axes[idx].scatter(x, y, alpha=0.5)
        axes[idx].plot(x, y_pred, color='blue', label='Best Fit', linestyle='--')

        # Add y=x line
        axes[idx].plot(x, x, color='black', linestyle='--', label='y = x')

        # Set titles and labels
        axes[idx].set_title(f'{cancer_type} (R² = {r_squared:.2f})')
        axes[idx].set_xlabel('-log10(p_val_HMF)')
        axes[idx].set_ylabel('-log10(p_val_PROFILE)')
        axes[idx].legend()
        axes[idx].set_xlim(0, max(x.max(), y.max()))
        axes[idx].set_ylim(0, max(x.max(), y.max()))

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

