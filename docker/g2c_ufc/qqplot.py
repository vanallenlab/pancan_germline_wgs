from qqman import qqman
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
	df_pancancer_assoc = pd.read_csv("/Users/noah/Desktop/UFC/saige_output/saige_gene_output.pancancer.txt",sep='\t')
	df_pancancer_assoc = df_pancancer_assoc[(df_pancancer_assoc['Group'] == 'Pathogenic;PLP;Likely_pathogenic;HIGH') & (df_pancancer_assoc['max_MAF'] == 0.01)]
	df_breast_assoc = pd.read_csv("/Users/noah/Desktop/UFC/saige_output/saige_gene_output.breast.txt",sep='\t')
	df_breast_assoc = df_breast_assoc[(df_breast_assoc['Group'] == 'Pathogenic;PLP;Likely_pathogenic;HIGH') & (df_breast_assoc['max_MAF'] == 0.01)]
	df_colorectal_assoc = pd.read_csv("/Users/noah/Desktop/UFC/saige_output/saige_gene_output.colorectal.txt",sep='\t')
	df_colorectal_assoc = df_colorectal_assoc[(df_colorectal_assoc['Group'] == 'Pathogenic;PLP;Likely_pathogenic;HIGH') & (df_colorectal_assoc['max_MAF'] == 0.01)]

	df_lung_assoc = pd.read_csv("/Users/noah/Desktop/UFC/saige_output/saige_gene_output.lung.txt",sep='\t')
	df_lung_assoc = df_lung_assoc[(df_lung_assoc['Group'] == 'Pathogenic;PLP;Likely_pathogenic;HIGH') & (df_lung_assoc['max_MAF'] == 0.01)]
	df_melanoma_assoc = pd.read_csv("/Users/noah/Desktop/UFC/saige_output/saige_gene_output.melanoma.txt",sep='\t')
	df_melanoma_assoc = df_melanoma_assoc[(df_melanoma_assoc['Group'] == 'Pathogenic;PLP;Likely_pathogenic;HIGH') & (df_melanoma_assoc['max_MAF'] == 0.01)]
	df_prostate_assoc = pd.read_csv("/Users/noah/Desktop/UFC/saige_output/saige_gene_output.prostate.txt",sep='\t')
	df_prostate_assoc = df_prostate_assoc[(df_prostate_assoc['Group'] == 'Pathogenic;PLP;Likely_pathogenic;HIGH') & (df_prostate_assoc['max_MAF'] == 0.01)]

	pancancer_p_vals = list(df_pancancer_assoc['Pvalue'])
	breast_p_vals = list(df_breast_assoc['Pvalue'])
	colorectal_p_vals = list(df_colorectal_assoc['Pvalue'])
	lung_p_vals = list(df_lung_assoc['Pvalue'])
	melanoma_p_vals = list(df_melanoma_assoc['Pvalue'])
	prostate_p_vals = list(df_prostate_assoc['Pvalue'])
	
	figure, axes = plt.subplots(nrows=2, ncols=3, figsize = (8,6))

	# qqman.qqplot("/Users/noah/Desktop/UFC/saige_output/saige_gene_output.breast.txt", ax=axes[0,0],title="From file")
	qqman.qqplot(pancancer_p_vals, ax=axes[0,0],title="Pancancer Cancer")
	qqman.qqplot(breast_p_vals, ax=axes[0,1],title="Breast Cancer")
	qqman.qqplot(colorectal_p_vals, ax=axes[0,2],title="Colorectal Cancer")
	qqman.qqplot(lung_p_vals, ax=axes[1,0],title="Lung Cancer")
	qqman.qqplot(melanoma_p_vals, ax=axes[1,1],title="Melanoma Cancer")
	qqman.qqplot(prostate_p_vals, ax=axes[1,2],title="Prostate Cancer")
	# qqman.qqplot(df_assoc.P, ax=axes[1,0],title="From Series")

	figure.tight_layout()
	plt.show()
	#plt.savefig("./SubQQplot.png",format="png")
	#plt.clf()
	#plt.close()