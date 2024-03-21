import os
from google.cloud import storage

# Initialize Google Cloud Storage client
storage_client = storage.Client()

# Specify bucket name and directory
bucket_name = 'fc-9cb68074-23c1-4bb3-9ef2-7363efd1fb40'
#directory = 'gsc'
directory = 'vaf_info'
# Initialize an empty string to store concatenated content
all_content = ''

# Filter blobs that match the pattern
blobs = storage_client.list_blobs(bucket_name, prefix=directory)
#filtered_blobs = [blob for blob in blobs if blob.name.endswith('.som_germ.txt') and blob.name.startswith(f'{directory}/HMF')]
filtered_blobs = [blob for blob in blobs if blob.name.endswith('.vaf.tsv') and blob.name.startswith(f'{directory}/HMF')]
# Iterate over each filtered blob
for blob in filtered_blobs:
    # Download blob content as string
    content = blob.download_as_text()

    # Append content to the 'all_content' string
    all_content += content

# Ensure the directory exists for the output file
output_directory = '/Users/noah/Desktop/DFCI_Data/gsi/data'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Write results to a new tab-delimited file
#output_file_path = os.path.join(output_directory, 'gsc.tsv')
output_file_path = os.path.join(output_directory, 'vaf_info.tsv')
with open(output_file_path, 'w') as f:
    f.write(all_content)

print(f'All files in {bucket_name}/{directory} matching the pattern are concatenated and saved to {output_file_path}.')
