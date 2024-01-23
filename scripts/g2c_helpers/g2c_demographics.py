import pandas as pd
import gcsfs
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import random

# File paths
local_tsv_path = '/Users/noah/Desktop/g2c_demo.tsv'

# Read the local TSV file
df = pd.read_csv(local_tsv_path, sep='\t')

cohorts = ['aou','biome','mesa','ufc','gtex','cptac','hgsvc','icgc','hmf','lcins','proactive-core','proactive-other']
ancestries = ['European', 'Latin-American-1','Latin-American-2','South-Asian','Asian-Pacific-Islander','East-Asian','African', 'African-American','Unknown','Other']
bar_colors = ['blue', 'red', 'pink', 'purple', 'darkgreen', 'lightgreen', 'yellow', 'orange', 'turquoise', 'brown']


# Dictionary to keep track of statistics
cohort_stats = {}
cohort_sex_stats = {}
overall_stats = {}
overall_sex_stats = {}
cohort_count = {}
total_count = 0

def grab_demographics_files():
    def list_demographics_files(base_path, cohort):
        # Create a GCS file system instance
        gcs_filesystem = gcsfs.GCSFileSystem(project='terra-446ad1f7')  # Replace with your actual project ID

        # Construct the path pattern for the specified cohort
        cohort_path = f"{base_path}/{cohort}/demographics/*"

        # List all files matching the pattern
        demographics_files = gcs_filesystem.glob(cohort_path)

        return demographics_files

    # Specify the base path in the "gs://dfci-g2c-inputs" file system and the cohort
    base_path = 'gs://dfci-g2c-inputs'
    
    demographics_files = []
    for cohort in cohorts:
        # Get a list of file pathways that have 'demographics' in them for the specified cohort
        demographics_files.extend([f"gs://{file}" for file in list_demographics_files(base_path, cohort)])

    return demographics_files

# Function to process each row
def process_file(demo_file):
    global total_count
    temp_df = pd.read_csv(demo_file,sep='\t')
    # Extract required columns
    temp_df = temp_df[['ancestry', 'chrX_count', 'chrY_count']]
    cohort = demo_file.split('/')[-3]
    cohort_count[cohort] = cohort_count.get(cohort, 0) + 1
    total_count += 1
    # Update cohort-specific statistics
    if cohort not in cohort_stats:
        cohort_stats[cohort] = {}
        cohort_sex_stats[cohort] = {}

    for _, temp_row in temp_df.iterrows():
        ancestry = temp_row['ancestry']
        chrX_count = int(temp_row['chrX_count'])
        chrY_count = int(temp_row['chrY_count'])
        sex = None
        if chrX_count == 1 and chrY_count == 1:
            sex = "Male"
        elif chrX_count == 2 and chrY_count == 0:
            sex = "Female"
        else:
            sex = "Aneuploidy"

        # Update cohort-specific statistics
        if ancestry not in cohort_stats[cohort]:
            cohort_stats[cohort][ancestry] = 0
        if sex not in cohort_sex_stats[cohort]:
            cohort_sex_stats[cohort][sex] = 0

        if ancestry not in overall_stats:
            overall_stats[ancestry] = 0
        if sex not in overall_sex_stats:
            overall_sex_stats[sex] = 0


        cohort_stats[cohort][ancestry] += 1
        cohort_sex_stats[cohort][sex] += 1
        overall_stats[ancestry] += 1
        overall_sex_stats[sex] += 1




# Loop through the rows of the DataFrame
for demo_file in tqdm(grab_demographics_files()):
    process_file(demo_file)


with open('g2c_demographics.txt', 'w') as file:
    # Write cohort-specific statistics
    file.write("=" * 30 + "\n")
    for cohort, stats in cohort_stats.items():
        file.write(f"Cohort: {cohort}\n")

        for ancestry in cohort_stats[cohort]:
            file.write(f"{ancestry}: {cohort_stats[cohort][ancestry]}\n")
        file.write("- " * 30 + "\n")
        for sex in cohort_sex_stats[cohort]:
            file.write(f"{sex}: {cohort_sex_stats[cohort][sex]}\n")
        file.write(f"Total: {cohort_count[cohort]}\n")
        file.write("=" * 30 + "\n")

    # Write overall statistics
    file.write("Overall Statistics\n")
    for ancestry in overall_stats:
        file.write(f"{ancestry}: {overall_stats[ancestry]}\n")
    file.write("- " * 30 + "\n")
    
    # Write overall sex statistics
    for sex in overall_sex_stats:
        file.write(f"{sex}: {overall_sex_stats[sex]}\n")

    file.write(f"Total: {total_count}\n")


def create_pie_chart(cohort,cohort_stats):
    values = list(cohort_stats[cohort].values())
    keys = list(cohort_stats[cohort].keys())

    #labels = [f"{key}: {value}" for key, value in zip(keys, values)]

    # Calculate the total count for percentage calculation
    total_count = sum(values)

    plt.pie(values,labels=keys,colors = bar_colors, autopct=lambda p: f'{p:.1f}%\n')
    plt.title(cohort, fontdict={'fontsize': 16, 'fontweight': 'bold'})
    plt.savefig(f"../../figures/g2c_figures/{cohort}_demographics.png")
    plt.clf()

def create_total_pie_chart(overall_stats):
    values = list(overall_stats.values())
    keys = list(overall_stats.keys())


    # Calculate the total count for percentage calculation
    total_count = sum(values)

    plt.pie(values,labels=keys,colors = bar_colors, autopct=lambda p: f'{p:.1f}%\n')
    plt.title("Everyone", fontdict={'fontsize': 16, 'fontweight': 'bold'})
    plt.savefig(f"../../figures/g2c_figures/Everyone_demographics.png")
    plt.clf()


# # Create subplots
# fig, ax = plt.subplots(figsize=(14, 8))

# Create bar graphs for each cohort and total
for cohort in cohort_stats:
    #print(cohort)
    if cohort not in cohort_stats.keys():
        continue
    #create_bar_graph(ax, cohort, cohort_stats, bar_width=0.2)
    create_pie_chart(cohort,cohort_stats)
create_total_pie_chart(overall_stats)

# # Add labels and title
# ax.set_xlabel('Ancestry Groups')
# ax.set_ylabel('Percentage')
# ax.set_title('Ancestry Distribution Across Cohorts')
# ax.set_xticks(np.arange(len(ancestries)) + 0.2)
# ax.set_xticklabels(ancestries)
# ax.legend()

# # Save the plot to 'g2c_demographics.png'
# plt.savefig('g2c_demographics.png')

# # Show the plot
# plt.show()




