#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


import argparse

def process_data(input_file, ancestry, ancestry_ID, chrX_count, chrY_count):
    # Read the file line by line
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Extract coordinates from the comments
    coordinates = [line.strip().split()[1:] for line in lines[3:6]]
    coordinate_dict = {item[0].rstrip(':'): tuple(map(float, item[1:])) for item in coordinates}


    # Extract data from the main part of the file
    data_line = lines[8].strip().split('\t')


    # Prepare output data
    output_data = {
        'Sample': data_line[0].strip(),
        'SNPs': int(data_line[1].strip()),
        'GD1(x)': float(data_line[2].strip()),
        'GD2(y)': float(data_line[3].strip()),
        'GD3(z)': float(data_line[4].strip()),
        'GD4': float(data_line[5].strip()),
        'E(%)': float(data_line[6].strip()),
        'F(%)': float(data_line[7].strip()),
        'A(%)': float(data_line[8].strip()),
        'F(xyz)': coordinate_dict['F'],
        'A(xyz)': coordinate_dict['A'],
        'E(xyz)': coordinate_dict['E'],
        'chrX_count': chrX_count,
        'chrY_count': chrY_count,
        'ancestry': ancestry,
        'ancestry_ID': ancestry_ID
    }

    # Convert output data to CSV format
    output_csv = (
        f"Sample\tSNPs\tGD1(x)\tGD2(y)\tGD3(z)\tGD4\tE(%)\tF(%)\tA(%)\tF(xyz)\tA(xyz)\tE(xyz)\tchrX_count\tchrY_count\tancestry\tancestry_ID\n"
        f"{output_data['Sample']}\t{output_data['SNPs']}\t{output_data['GD1(x)']}\t{output_data['GD2(y)']}\t"
        f"{output_data['GD3(z)']}\t{output_data['GD4']}\t{output_data['E(%)']}\t{output_data['F(%)']}\t{output_data['A(%)']}\t"
        f"{output_data['F(xyz)']}\t{output_data['A(xyz)']}\t{output_data['E(xyz)']}\t{output_data['chrX_count']}\t"
        f"{output_data['chrY_count']}\t{output_data['ancestry']}\t{output_data['ancestry_ID']}\n"
    )

    print(output_csv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process input data and create output TSV.')
    parser.add_argument('--file', type=str, help='Input file path', required=True)
    parser.add_argument('--ancestry', type=str, help='Column name for Ancestry', required=True)
    parser.add_argument('--ancestryID', type=int, help='Column name for Ancestry ID', required=True)
    parser.add_argument('--chrX_count', type=int, help='Column name for chrX count', required=True)
    parser.add_argument('--chrY_count', type=int, help='Column name for chrY count', required=True)

    args = parser.parse_args()

    process_data(args.file, args.ancestry, args.ancestryID, args.chrX_count, args.chrY_count)
