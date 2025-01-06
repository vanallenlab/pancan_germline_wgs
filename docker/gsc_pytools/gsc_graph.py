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
        print(f"{cancer_type}: has {len(data)} rows")


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
                 germline_event="germline_risk_allele",
                 figure_title="Volcano Plot of Germline Somatic Convergences",
                 save_path="insert_path_here.png"):

    # Drop rows where either p-value or OR is NaN
    df_clean = df.dropna(subset=[p_col, or_col])
    df = df_clean
    
    # Calculate -log10(p-value) and log2(OR)
    df['log10_p'] = -np.log10(df[p_col])
    df['log2_OR'] = np.log2(df[or_col])
    
    # Define shapes and colors
    shapes = {'direct_ppi_known_function': 'o',
             'direct_ppi_known_function;direct_ppi_unspecified_function':'o',
             'direct_ppi_known_function;nonspecific_ppi':'o',
             'direct_ppi_known_function;nonspecific_ppi;protein_complex':'o',
             'direct_ppi_unspecified_function':'h',
             'direct_ppi_unspecified_function;nonspecific_ppi':'h',
             'direct_ppi_unspecified_function;ligand_receptor':'^',
             'direct_ppi_unspecified_function;ligand_receptor;nonspecific_ppi':'^',
             'direct_ppi_unspecified_function;ligand_receptor;nonspecific_ppi;protein_complex':'^',
             'ligand_receptor':'^',
             'ligand_receptor;nonspecific_ppi':"^",
             'protein_complex':'h',
             'nonspecific_ppi;protein_complex':'h',
             'direct_ppi_unspecified_function;protein_complex':'h',
             'direct_ppi_unspecified_function;nonspecific_ppi;protein_complex':'h',
             'same_gene': 's',
             'nonspecific_ppi':'*'}
    colors = {'Breast': 'red', 'Prostate': 'blue', 'Colorectal': 'green', 
              'Lung': 'purple', 'Kidney': 'orange', 'Pancancer': 'black'}
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot points
    for criteria_val, shape in shapes.items():
        for cancer_type, color in colors.items():
            subset = df[(df[criteria_col] == criteria_val) & 
                        (df[cancer_type_col] == cancer_type) & 
                        (df['log2_OR'] > -20)]
            ax.scatter(subset['log2_OR'], subset['log10_p'], 
                       marker=shape, color=color, alpha=0.7)

    # Add significance thresholds
    bonferroni_threshold = -np.log10(0.05 / len(df))
    nominal_threshold = -np.log10(0.05)
    rejected, pvals_corrected = fdrcorrection(df[p_col], alpha=0.05)
    if rejected.any():
        fdr_threshold = -np.log10(pvals_corrected[rejected].max())
        ax.axhline(y=fdr_threshold, color='green', linestyle='--', label=f'FDR Threshold')
    ax.axhline(y=bonferroni_threshold, color='red', linestyle='--', label=f'Bonferroni')
    ax.axhline(y=nominal_threshold, color='blue', linestyle='--', label='Nominal')

    # Annotate points above Bonferroni
    significant_points = df[df['log10_p'] > nominal_threshold]
    for _, row in significant_points.iterrows():
        # Get the value from germline_event, fallback to 'germline_gene' if germline_event is NA or NaN
        germline_value = row[germline_event] if pd.notna(row[germline_event]) else row['germline_gene']
        label_text = f"({germline_value}, {row['somatic_gene']})"
        ax.text(row['log2_OR'], row['log10_p'], label_text, fontsize=8, ha='right')

    # # Add legends
    # color_handles = [plt.Line2D([0], [0], marker='o', color='w', 
    #                             markerfacecolor=color, label=cancer) 
    #                  for cancer, color in colors.items()]
    # color_legend = ax.legend(handles=color_handles, title="Cancer Type", 
    #                          bbox_to_anchor=(1.05, 1), loc='upper left')
    # ax.add_artist(color_legend)

    # shape_handles = [plt.Line2D([0], [0], marker=shape, color='k', label=criteria) 
    #                  for criteria, shape in shapes.items()]
    # shape_legend = ax.legend(handles=shape_handles, title="Criteria", 
    #                          bbox_to_anchor=(1.05, 0.75), loc='upper left')
    # ax.add_artist(shape_legend)

    # Add vertical dashed line at x=0
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, label='x=0')

    # # # Dashed lines legend
    # # line_handles = [
    # #     plt.Line2D([0], [0], color='red', linestyle='--', label='Bonferroni'),
    # #     plt.Line2D([0], [0], color='blue', linestyle='--', label='Nominal'),
    # #     plt.Line2D([0], [0], color='green', linestyle='--', label='Benjamini-Hochberg'),
    # #     plt.Line2D([0], [0], color='gray', linestyle='--', label='x=0')
    # # ]

    # # # Use a separate legend for lines
    # # line_legend = ax.legend(handles=line_handles, title="Significance Lines", 
    # #                     bbox_to_anchor=(1.05, 0.5), loc='upper left')
    # # ax.add_artist(line_legend)

    # Add point legends (color and shape combined)
    color_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, label=cancer) 
        for cancer, color in colors.items()
    ]
    shape_handles = [
        plt.Line2D([0], [0], marker=shape, color='k', label=criteria) 
        for criteria, shape in shapes.items()
    ]
    
    # Add line legend for thresholds
    line_handles = [
        plt.Line2D([0], [0], color='red', linestyle='--', label='Bonferroni'),
        plt.Line2D([0], [0], color='blue', linestyle='--', label='Nominal'),
        plt.Line2D([0], [0], color='green', linestyle='--', label='FDR Threshold'),
        plt.Line2D([0], [0], color='gray', linestyle='--', label='x=0'),
    ]

    # Create a combined legend outside the plot
    fig.legend(handles=color_handles + shape_handles + line_handles,
               title="Legend", loc='center right', bbox_to_anchor=(1.3, 0.5))

    # Adjust layout to give room for the legend
    plt.subplots_adjust(right=0.75)
    
    # Labels, title, and save
    ax.set_xlabel('log2(OR)')
    ax.set_ylabel('-log10(p)')
    ax.set_title(figure_title)
    plt.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)



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

def plot_gwas_subplots(df, save_path='gwas_subplots.png'):
    df = df[(df['relevant_cancer'] == 1) | (df['relevant_cancer'] == 2)]
    # Convert p-values to -log10 scale for plotting
    df['-log10_p_val_HMF'] = -np.log10(df['p_val_HMF'])
    df['log2OR'] = np.log2(df['OR_HMF'])
    # Extract chromosome position as integer from `germline_risk_allele`
    df['chrom_position'] = df['germline_risk_allele'].apply(lambda x: int(x.split(':')[1].split('-')[0]))

    # Group by cancer_type, germline_gene, and somatic_gene
    grouped = df.groupby(['cancer_type', 'germline_gene', 'somatic_gene'])

    # Calculate nominal and Bonferroni thresholds
    nominal_threshold = -np.log10(0.05)

    # Set up the plot grid
    num_plots = len(grouped)
    fig, axes = plt.subplots(nrows=(num_plots // 3) + 1, ncols=3, figsize=(15, 5 * (num_plots // 3 + 1)))
    axes = axes.flatten()

    for i, (name, group) in enumerate(grouped):
        cancer_type, germline_gene, somatic_gene = name
        print(cancer_type)
        # Calculate Bonferroni threshold for the current group
        bonferroni_threshold = -np.log10(0.05 / len(group))

        # Plot on the i-th subplot
        ax = axes[i]
        ax.scatter(group['chrom_position'], group['-log10_p_val_HMF'], color='b', s=10, alpha=0.6)
        #ax.scatter(group['chrom_position'], group['log2OR'], color='b', s=10, alpha=0.6)
        ax.axhline(nominal_threshold, color='green', linestyle='--', label='Nominal p=0.05')
        ax.axhline(bonferroni_threshold, color='red', linestyle='-', label='Bonferroni Corrected')

        # Set plot titles and labels
        title = ""
        if cancer_type == "Pancancer":
            title = f"{germline_gene} Noncoding and {somatic_gene} Coding Convergence: Pancancer"
        else:
            title = f"{germline_gene} Noncoding and {somatic_gene} Coding Convergence: {cancer_type}"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Chromosome Position")
        ax.set_ylabel("-log10(p)")
        ax.legend()

    # Remove any extra subplots in the grid
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.show()
    print(f"Figure saved to {save_path}")