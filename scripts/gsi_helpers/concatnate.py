#Concatnate all files from one gcs bucket into one -input parameters needed.
import os
from google.cloud import storage

# Initialize Google Cloud Storage client
storage_client = storage.Client()

# Specify bucket name and directory
bucket_name = 'fc-9cb68074-23c1-4bb3-9ef2-7363efd1fb40'
directory = 'gsc'

# Initialize an empty string to store concatenated content
all_content = ''

# Iterate over each blob in the specified directory
blobs = storage_client.list_blobs(bucket_name, prefix=directory)
for blob in blobs:
    # Download blob content as string
    content = blob.download_as_text()
    
    # Append content to the 'all_content' string
    all_content += content

# Upload the concatenated content to a new blob
# destination_blob_name = f'{directory}/hiv_all_samples.txt'
# bucket = storage_client.bucket(bucket_name)
# blob = bucket.blob(destination_blob_name)
# blob.upload_from_string(all_content)

# Ensure the directory exists for the output file
output_directory = 'data'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Write results to a new tab-delimited file
output_file_path = os.path.join(output_directory, 'gsc.tsv')
with open(output_file_path, 'w') as f:
    f.write(all_content)

print(f'All files in {bucket_name}/{directory} are concatenated and saved to {output_file_path}.')
