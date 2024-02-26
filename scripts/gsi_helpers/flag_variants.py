from google.cloud import storage
import os

# Initialize Google Cloud Storage client
storage_client = storage.Client()

# Specify bucket name and file path
bucket_name = 'fc-9cb68074-23c1-4bb3-9ef2-7363efd1fb40'
file_path = 'hiv_tracker/hiv_all_samples.txt'

# Function to process the content of hiv_all_samples.txt
def flag_variants(content):
    # Count occurrences of each value
    value_counts = {}
    #lines = set()
    for line in content.split('\n'):
        #lines.add(line)
        try:
            values = line.split('\t')[2] + ":" + line.split('\t')[3]
        except:
            print(line)


        value_counts[values] = value_counts.get(values, 0) + 1
    # Calculate thresholds
    flag_threshold = 0.01 * 2115
    warning_threshold = 0.05 * 2115
    failure_threshold = 0.10 * 2115
    
    # Determine values exceeding thresholds
    results = []
    for val, count in value_counts.items():
        if count > failure_threshold:
            results.append((val, 'failure',count))
        elif count > warning_threshold:
            results.append((val, 'warning',count))
        elif count > flag_threshold:
            results.append((val,'flag',count))
        else:
            results.append((val,'good',count))

    return results



# Download hiv_all_samples.txt from Google Cloud Storage
# bucket = storage_client.bucket(bucket_name)
# blob = bucket.blob(file_path)
# content = blob.download_as_text()

# Open the file in read mode ('r')
with open('data/gsc.tsv', 'r') as file:
    # Read the entire contents of the file as one long string
    content = file.read()
# Process the content
results = flag_variants(content)


# Ensure the directory exists for the output file
output_directory = 'data'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Write results to a new tab-delimited file
output_file_path = os.path.join(output_directory, 'warnings_and_failures.tsv')
with open(output_file_path, 'w') as f:
    # Write header
    f.write("gene_loc\tflag\tcount\n")
    # Write results
    for val, tag, count in results:
        f.write(f"{val}\t{tag}\t{count}\n")


print(f'Results written to {output_file_path}.')
