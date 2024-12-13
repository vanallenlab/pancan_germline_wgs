from qqman import qqman
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    df_pancancer_assoc = pd.read_csv("/Users/noah/Desktop/UFC/saige_output/saige_gene_output.pancancer.txt", sep='\t')
    df_pancancer_assoc1 = df_pancancer_assoc[
        (df_pancancer_assoc['Group'] == 'Pathogenic;PLP;Likely_pathogenic;HIGH') &
        (df_pancancer_assoc['max_MAF'] == 0.01)
    ]
    df_pancancer_assoc2 = df_pancancer_assoc[
        (df_pancancer_assoc['Group'] == 'Pathogenic;PLP;Likely_pathogenic;HIGH;MODERATE') &
        (df_pancancer_assoc['max_MAF'] == 0.01)
    ]

    pancancer_p_vals = list(df_pancancer_assoc['Pvalue'])
    
    figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(30, 15))

    # First plot: High Impact Variants
    qqman.qqplot(df_pancancer_assoc1.Pvalue, ax=axes[0])
    axes[0].set_title("PLP and High Impact Variants", fontsize=44, pad=30)
    axes[0].tick_params(axis='both', labelsize=40)  # Adjust tick size
    axes[0].set_xlabel("Expected -log10(p)", fontsize=40, weight='bold')
    axes[0].set_ylabel("Observed -log10(p)", fontsize=40, weight='bold')

    # Second plot: High/Moderate Impact Variants
    qqman.qqplot(df_pancancer_assoc2.Pvalue, ax=axes[1])
    axes[1].set_title("PLP and High/Moderate Impact Variants", fontsize=44, pad=30)
    axes[1].tick_params(axis='both', labelsize=40)  # Adjust tick size
    axes[1].set_xlabel("Expected -log10(p)", fontsize=40, weight='bold')
    axes[1].set_ylabel("Observed -log10(p)", fontsize=40, weight='bold')

    # Add overarching title
    figure.suptitle("Preliminary Gene-Based Rare Variant Association Testing on Chr22", fontsize=50, weight='bold')  # Adjust y to bring title down

    # Add footnote with more space (lower the y position)
    figure.text(0.5, 0.01, "*PLP = Pathogenic and/or Likely Pathogenic", ha='center', fontsize=30, style='italic')

    # Adjust spacing between plots
    plt.subplots_adjust(wspace=0.8, top=0.83, bottom=0.2)  # Adjust horizontal space between subplots

    # Save the figure
    plt.savefig("/Users/noah/Desktop/UFC/plots/QQplot_BroadRetreat_2024.png", format="png", dpi=800)
