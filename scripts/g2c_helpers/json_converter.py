# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

import json
import argparse
import os
import re

def process_json(input_file, output_file):
    line = open(input_file).readline()[1:-1]
    line = line.replace('\":\"','\" : \"')
    line = re.sub(r'\${(.*?)}', r'\1', line)
    line = line.replace('\"true\"','true')
    line = line.replace('\"false\"','false')
    line = line.replace('\"[','[').replace(']\"',']')

    arr = re.split(r',(?![^\[]*\])', line)

    with open(output_file, 'w') as f:
        f.write('{\n')
        for val in arr:
            val = re.sub(r'"(\d+(\.\d+)?)",?', r'\1', val)
            if val != arr[-1]:
                f.write(f"\t{val},\n")
            else:
                f.write(f"\t{val}\n")
        f.write('}')
    print(f"Output written to: {output_file}")



def main():
    parser = argparse.ArgumentParser(description='Convert JSON file keys to separate lines.')
    parser.add_argument('input_file', help='Input JSON file')
    parser.add_argument('output_file', help='Output JSON file')

    args = parser.parse_args()

        
    process_json(args.input_file, args.output_file)

if __name__ == '__main__':
    main()

