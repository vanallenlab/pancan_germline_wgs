import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

def plot_germline_frequencies(data, 
                              hmf_col='germline_plp_frequency_HMF', 
                              profile_col='germline_plp_frequency_PROFILE', 
                              save_path="germline_snp_frequencies.png"):

    # Define the cancer types of interest
    cancer_types = ['Breast', 'Prostate', 'Colorectal', 'Lung', 'Kidney', 'Pancancer']
    
    # Drop rows where either the HMF or PROFILE frequency is missing (NaN)
    data_clean = data.dropna(subset=[hmf_col, profile_col])
    data = data_clean

    # Create a figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()  # Flatten the array for easy indexing

    # Iterate through each cancer type and create a subplot
    for idx, cancer_type in enumerate(cancer_types):

        # Filter the DataFrame for the current cancer type
        filtered_data = data[data['cancer_type'] == cancer_type]

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
        axes[idx].set_xlabel('germline_plp_frequency_HMF')
        axes[idx].set_ylabel('germline_plp_frequency_PROFILE')
        axes[idx].legend()
        axes[idx].set_xlim(0, max(x.max(), y.max()))
        axes[idx].set_ylim(0, max(x.max(), y.max()))

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_somatic_mutation_frequencies(df, 
                                        x_col="somatic_mutation_frequency_HMF", 
                                        y_col="somatic_mutation_frequency_PROFILE", 
                                        save_path="somatic_mutation_frequencies.png"):
    # Specify the cancer types to include
    cancer_types = ["Breast", "Prostate", "Colorectal", "Lung", "Kidney", "Pancancer"]

    # Drop rows where either the HMF or PROFILE frequency is missing (NaN)
    df_clean = df.dropna(subset=[x_col, y_col])
    df = df_clean

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
                 save_path="volcano_plot.png"):
    

    # Drop rows where either the HMF or PROFILE frequency is missing (NaN)
    df_clean = df.dropna(subset=[p_col, or_col])
    df = df_clean
    
    # Calculate -log10(p-value) and log2(OR)
    df['log10_p'] = -np.log10(df[p_col])
    df['log2_OR'] = np.log2(df[or_col])
    
    # Set up point shapes based on 'criteria' and colors based on 'cancer_type'
    shapes = {'known_ppi': 'o', 'other': 's'}  # Example shapes for criteria
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
    plt.axhline(y=bonferroni_threshold, color='red', linestyle='--', label='Bonferroni Significance')
    plt.axhline(y=nominal_threshold, color='blue', linestyle='--', label='Nominal Significance')
    
    # Label points above Bonferroni threshold with one of the 5 lowest p-values
    significant_points = df[df['log10_p'] > bonferroni_threshold].nsmallest(5, p_col)
    for _, row in significant_points.iterrows():
        plt.text(row['log2_OR'], row['log10_p'], row[or_col], fontsize=8, ha='right')

    # Set plot labels and title
    plt.xlabel('log2(OR_combined)')
    plt.ylabel('-log10(p_val_combined)')
    plt.title('Volcano Plot with Bonferroni and Nominal Significance')
    
    # Display legend
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Save plot to file
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
