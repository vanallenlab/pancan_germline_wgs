import pandas as pd

# Define file paths
control_file = 'superset_control_demographics.tsv'
case_file = 'case_demographics.tsv'
output_file = 'control_demographics_interim.tsv'

# Read case and control files
control_df = pd.read_csv(control_file, sep='\t', index_col=0)
case_df = pd.read_csv(case_file, sep='\t', index_col=0)

# Sort control and case data
sorted_control = control_df.sort_values(by=['race', 'ethnicity', 'year_of_birth', 'sex_at_birth'])
sorted_case = case_df.sort_values(by=['race', 'ethnicity', 'year_of_birth', 'sex_at_birth'])

# Initialize matched controls dictionary
matched_controls = {}
print(sorted_control.shape,sorted_case.shape)
# Exact Matching
for idx, case_row in sorted_case.iterrows():
    match = sorted_control[
        (sorted_control['race'] == case_row['race']) &
        (sorted_control['ethnicity'] == case_row['ethnicity']) &
        (sorted_control['year_of_birth'] == case_row['year_of_birth']) &
        (sorted_control['sex_at_birth'] == case_row['sex_at_birth'])
    ].head(1)
    if not match.empty:        control_idx = match.index[0]        matched_controls[idx] = control_idx #Case is key and control is value        sorted_control = sorted_control.drop(control_idx) #Drops control value from controls df
        sorted_case = sorted_case.drop(idx)
#This is all four variables of interest.
for bandwidth in range(1,20):
    for idx, case_row in sorted_case.iterrows():
        #print(f"Bandwidth: {bandwidth}")
        match = sorted_control[
            (sorted_control['race'] == case_row['race']) &
            (sorted_control['ethnicity'] == case_row['ethnicity']) &
            (sorted_control['year_of_birth'].between(case_row['year_of_birth'] - bandwidth, case_row['year_of_birth'] + bandwidth)) &
            (sorted_control['sex_at_birth'] == case_row['sex_at_birth'])
        ].head(1)
        if not match.empty:
            control_idx = match.index[0]
            matched_controls[idx] = control_idx #Case is key and control is value            sorted_control = sorted_control.drop(control_idx) #Drops control value from controls df            sorted_case = sorted_case.drop(idx)
#This is just for race, age, and sex. We exclude ethnicity because that is often tied to race.
for bandwidth in range(1,20):
    for idx, case_row in sorted_case.iterrows():
        #print(f"Bandwidth: {bandwidth}")
        match = sorted_control[
            (sorted_control['race'] == case_row['race']) &
            (sorted_control['year_of_birth'].between(case_row['year_of_birth'] - bandwidth, case_row['year_of_birth'] + bandwidth)) &
            (sorted_control['sex_at_birth'] == case_row['sex_at_birth'])
        ].head(1)
        if not match.empty:
            control_idx = match.index[0]
            matched_controls[idx] = control_idx #Case is key and control is value
            sorted_control = sorted_control.drop(control_idx) #Drops control value from controls df
            sorted_case = sorted_case.drop(idx)

#This is just for race and age. We are putting race ahead of sex for matching.
for bandwidth in range(1,20):
    for idx, case_row in sorted_case.iterrows():
        match = sorted_control[
            (sorted_control['race'] == case_row['race']) &
            (sorted_control['year_of_birth'].between(case_row['year_of_birth'] - bandwidth, case_row['year_of_birth'] + bandwidth))
        ].head(1)
        if not match.empty:
            control_idx = match.index[0]
            matched_controls[idx] = control_idx #Case is key and control is value
            sorted_control = sorted_control.drop(control_idx) #Drops control value from controls df
            sorted_case = sorted_case.drop(idx)

#This is just for sex and age once race and ethnicity have been used up
for bandwidth in range(1,20):
    for idx, case_row in sorted_case.iterrows():
        match = sorted_control[
            (sorted_control['sex_at_birth'] == case_row['sex_at_birth']) &
            (sorted_control['year_of_birth'].between(case_row['year_of_birth'] - bandwidth, case_row['year_of_birth'] + bandwidth))
        ].head(1)
        if not match.empty:
            control_idx = match.index[0]
            matched_controls[idx] = control_idx #Case is key and control is value
            sorted_control = sorted_control.drop(control_idx) #Drops control value from controls df
            sorted_case = sorted_case.drop(idx)

#We are only looking for year of birth as a last resort
#This is useful for when all of the other fields are missing
for bandwidth in range(1,20):
    for idx, case_row in sorted_case.iterrows():
        match = sorted_control[(sorted_control['year_of_birth'].between(case_row['year_of_birth'] - bandwidth, case_row['year_of_birth'] +bandwidth))
        ].head(1)
        if not match.empty:
            control_idx = match.index[0]
            matched_controls[idx] = control_idx #Case is key and control is value
            sorted_control = sorted_control.drop(control_idx) #Drops control value from controls df
            sorted_case = sorted_case.drop(idx)





sorted_control.to_csv('unmatched.tsv',sep='\t')
sorted_case.to_csv('unmatchables.tsv',sep='\t')
print(len(sorted_case))

# Convert dictionary to DataFrame
df = pd.DataFrame(list(matched_controls.items()), columns=['person_id_case', 'person_id_control'])

# Save DataFrame to TSV file
output_file = 'matched_controls.tsv'
df.to_csv(output_file, sep='\t', index=False)

print("Exact matching done. Interim results saved to:", output_file)
