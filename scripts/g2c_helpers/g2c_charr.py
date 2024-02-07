# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

import pandas as pd
import gcsfs
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import random


cohorts = ['aou','biome','mesa','ufc','gtex','cptac','hgsvc','icgc','hmf','lcins','proactive-core','proactive-other']
ancestries = ['European', 'Latin-American-1','Latin-American-2','South-Asian','Asian-Pacific-Islander','East-Asian','African', 'African-American','Unknown','Other']
bar_colors = ['blue', 'red', 'pink', 'purple', 'darkgreen', 'lightgreen', 'yellow', 'orange', 'turquoise', 'brown']


# Dictionary to keep track of statistics
charr_arr = []
charr_reblocked_arr = []
charr_dict = dict()
variable_of_interest = 'CHARR'

def grab_charr_files():
    def list_charr_files(base_path, cohort):
        # Create a GCS file system instance
        gcs_filesystem = gcsfs.GCSFileSystem(project='terra-446ad1f7')  # Replace with your actual project ID

        # Construct the path pattern for the specified cohort
        cohort_path = f"{base_path}/{cohort}/charr/*"

        # List all files matching the pattern
        charr_files = gcs_filesystem.glob(cohort_path)

        return charr_files

    # Specify the base path in the "gs://dfci-g2c-inputs" file system and the cohort
    base_path = 'gs://dfci-g2c-inputs'
    
    charr_files = []
    for cohort in cohorts:
        # Get a list of file pathways that have 'demographics' in them for the specified cohort
        charr_files.extend([f"gs://{file}" for file in list_charr_files(base_path, cohort)])

    return charr_files

def grab_charr_reblocked_files():
    path = '/Users/noah/Desktop/ssi.tsv'
    df = pd.read_csv(path,sep='\t')
    charr_reblocked_files = list(df['contamination_reblocked'])

    return charr_reblocked_files

# Function to process each row
def process_file(charr_file):
    df = pd.read_csv(charr_file,sep='\t')
    return float(df[variable_of_interest].iloc[0])

def process_cohorts():
    charr_files = grab_charr_files()
    cohort_dict = dict()
    cohort_dict['mesa'] = [file for file in charr_files if 'mesa' in file]
    cohort_dict['gtex'] = [file for file in charr_files if 'gtex' in file]
    cohort_dict['icgc'] = [file for file in charr_files if 'icgc' in file]
    cohort_dict['lcins'] = [file for file in charr_files if 'lcins' in file]
    #cohort_dict['hmf'] = [file for file in charr_files if 'hmf' in file]
    cohort_dict['ufc'] = [file for file in charr_files if 'ufc' in file]
    cohort_dict['cptac'] = [file for file in charr_files if 'cptac' in file]
    cohort_dict['hgsvc'] = [file for file in charr_files if 'hgsvc' in file]
    cohort_dict['biome'] = [file for file in charr_files if 'biome' in file]
    cohort_dict['proactive-core'] = [file for file in charr_files if 'proactive-core' in file]
    cohort_dict['proactive-other'] = [file for file in charr_files if 'proactive-other' in file]


    
    for key in cohort_dict:
        for file in cohort_dict[key]:
            if key not in charr_dict:
                charr_dict[key]=[]
            charr_dict[key].append(process_file(file))

def make_figure():
    def box_plot(data,position, label):
        plt.boxplot(data,positions=[position], labels=[label])


    box_plot(charr_dict['mesa'],1,"MESA")
    box_plot(charr_dict['gtex'],2,"GTEX")
    box_plot(charr_dict['icgc'],3,"ICGC")
    box_plot(charr_dict['lcins'],4,"LCINS")
    #box_plot(charr_dict['hmf'],5,"HMF")
    box_plot(charr_dict['ufc'],6,"UFC")
    box_plot(charr_dict['cptac'],7,"CPTAC")
    box_plot(charr_dict['hgsvc'],8,"HGSVC")
    box_plot(charr_dict['biome'],9,"BIOME")
    box_plot(charr_dict['proactive-core'],10,"P-CORE")
    box_plot(charr_dict['proactive-other'],11,"P-OTHER")

    # Show the plot
    # Add labels and title
    plt.xlabel('Cohort')
    plt.ylabel('CHARR')
    plt.title(f'Box Plot of {variable_of_interest} values across Cohorts')

    # Add dotted lines and labels
    if variable_of_interest == "CHARR":
        plt.axhline(y=0.02, linestyle='--', color='gray', label='Subtle Contamination')
        plt.axhline(y=0.03, linestyle='--', color='red', label='Suggestive of Contamination')
        plt.legend()
    # Create a wider figure
    #plt.figure(figsize=(10, 6))  # Adjust width as needed
    plt.savefig(f"../../figures/charr_figures/{variable_of_interest}.png")

def compare_reblocked_raw():
    # Function to process each row
    def process_file_compare(charr_file):
        df = pd.read_csv(charr_file,sep='\t')
        return float(df['CHARR'].iloc[0])
    def list_charr_files(base_path):
        # Create a GCS file system instance
        gcs_filesystem = gcsfs.GCSFileSystem(project='terra-446ad1f7')  # Replace with your actual project ID

        # Construct the path pattern for the specified cohort
        original_path = f"{base_path}/charr_original/*"
        reblocked_path = f"{base_path}/charr/*"

        # List all files matching the pattern
        original_files = sorted(gcs_filesystem.glob(original_path))
        reblocked_files = [file.replace('charr_original','charr') for file in original_files]

        return original_files,reblocked_files

    # Specify the base path in the "gs://dfci-g2c-inputs" file system and the cohort
    base_path = 'gs://dfci-g2c-inputs/mesa'
    original_files,reblocked_files = list_charr_files(base_path)

    coordinates = []
    for i in range(len(original_files)):
        coordinates.append((process_file_compare("gs://" + original_files[i]),process_file_compare("gs://" + reblocked_files[i])))

    # Extract the first and second integers from each tuple into separate arrays
    x = np.array([t[0] for t in coordinates])
    y = np.array([t[1] for t in coordinates])

    correlation_coefficient = np.corrcoef(x, y)[0,1]
    r_squared = correlation_coefficient ** 2
    print(r_squared)



if __name__ == '__main__':
    compare_reblocked_raw()
    #process_cohorts()
    #make_figure()
    # main()


