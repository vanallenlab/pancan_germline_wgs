# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# The purpose of this script is to convert between json
# code from terra, put it into human readable form (that can be edited if necessary),
# and then convert to a json file that is ideal for AllofUs workbench.
# If you are starting from the original terra json, use the flag '-t',
# otherwise run w/ no flags.
import json
import argparse
import os
import re

# The goal of this function is to turn terra json
# files into a human readable json w/ the suffix .hread.json
# that can be easily edited.
def terra_to_readable_json(input_file, output_file):
    line = open(input_file).readline()[1:-1]
    line = line.replace('\":\"','\" : \"')
    line = re.sub(r'\${(.*?)}', r'\1', line)
    line = line.replace('\"true\"','true')
    line = line.replace('\"false\"','false')
    line = line.replace('\"[','[').replace(']\"',']')

    arr = re.split(r',(?![^\[]*\])', line)
    arr = [s for s in arr if '\"\"' not in s]
    arr = [s for s in arr if '[]' not in s]
    print(len(arr))

    with open(output_file, 'w') as f:
        f.write('{\n')
        for val in arr:

            val = re.sub(r'"(\d+(\.\d+)?)",?', r'\1', val) #turns strings into floats or ints if applicable
            val = val.replace('\\',"")
            if val != arr[-1]:
                f.write(f"\t{val},\n")
            else:
                f.write(f"\t{val}\n")
                
        f.write('}')
    
    print(f"Output written to: {output_file}")

# This function takes in an input file and
# writes to the directory of the output_file.
# This takes a human readable json file (normal json)
# and converts to a json file used by AllofUs workbench.
def readable_to_aou_json(input_file, output_file):
    try:
        # Open the input file for reading
        with open(input_file, 'r') as input_file_handle:
            # Read the content from the input file
            content = input_file_handle.read()

            # Remove '\n' and '\t' characters from the content
            content_without_newline_tabs = content.replace('\n', '').replace('\t', '')

            # Open the output file for writing
            with open(output_file, 'w') as output_file_handle:
                # Write the modified content to the output file
                output_file_handle.write(content_without_newline_tabs)

        print(f"Conversion successful. Output written to {output_file}")

    except Exception as e:
        print(f"Error: {e}")

# Example usage:
# readable_to_aou_json("input.txt", "output.txt")





def main():
    parser = argparse.ArgumentParser(description='Convert JSON file keys to separate lines.')
    parser.add_argument('input_file', help='Input JSON file')
    parser.add_argument('-t', '--terra', action='store_true', help='If present, run both functions')

    args = parser.parse_args()

    if args.terra:
        terra_to_readable_json(args.input_file, args.input_file.replace('.json','.hread.json'))
        readable_to_aou_json(args.input_file.replace('.json','.hread.json'), args.input_file.replace('.json','.out.json'))
    else:
        if '.hread.json' in args.input_file:
            readable_to_aou_json(args.input_file,args.input_file.replace('.hread.json','.out.json'))
        else:
            readable_to_aou_json(args.input_file,args.input_file.replace('.json','.out.json'))

if __name__ == '__main__':
    main()

