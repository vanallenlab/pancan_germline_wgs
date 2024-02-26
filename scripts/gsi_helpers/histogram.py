import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/hivs_germline_genes.tsv',sep='\t',names=['gene','label','count'])

plt.hist([df['count']])
plt.xlabel("Mutated in Cohort")
plt.ylabel("Frequency")
plt.show()