#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

import sys
import csv

def read_sample_data(filename):
    # Dictionary to store the data
    sample_data = {}

    # Read the text file and populate the dictionary
    with open(filename, 'r') as txt_file:
        lines = txt_file.readlines()

        # Skip lines until the header line
        header_index = lines.index('Sample\t#SNPs\tGD1 (x)\tGD2 (y)\tGD3 (z)\tGD4\tE(%)\tF(%)\tA(%)\n')

        for line in lines[header_index + 1:]:
            # Split the line into columns
            columns = line.strip().split('\t')

            # Extract relevant data
            sample = columns[0]
            gd1 = float(columns[2])
            gd2 = float(columns[3])
            gd3 = float(columns[4])
            gd4 = float(columns[5])
            e_percent = float(columns[6])
            f_percent = float(columns[7])
            a_percent = float(columns[8])

            # Store the values in the dictionary
            sample_data[sample] = {
                'GD1': gd1,
                'GD2': gd2,
                'GD3': gd3,
                'GD4': gd4,
                'E(%)': e_percent,
                'F(%)': f_percent,
                'A(%)': a_percent
            }

    return sample_data

def classify_population(sample_data):
    for sample, data in sample_data.items():
        e_percent = data['E(%)']
        f_percent = data['F(%)']
        a_percent = data['A(%)']
        gd1 = data['GD1']

        if e_percent >= 90:
            print(f"European\n1")
        elif f_percent >= 95:
            print(f"African\n2")
        elif a_percent >= 95:
            print(f"East-Asian\n3")
        elif 40 <= f_percent < 95 and a_percent < 14:
            print(f"African-American\n4")
        elif f_percent < 50 and e_percent < 90 and a_percent < 14 and gd1 < 1.48:
            print(f"Latin-American-1\n5")
        elif a_percent >= 14 and f_percent >= 14:
            print(f"Other\n9")
        elif gd1 > 30 * data['GD4']**2 + 1.58 and f_percent < 14:
            print(f"Asian-Pacific-Islander\n7")
        elif gd1 + data['GD4'] < 1.525 and f_percent < 14:
            print(f"Latin-American-2\n6")
        elif data['GD4'] > 5 * (gd1 - 1.524)**2 + 0.0575:
            print(f"South-Asian\n8")
        else:
            print(f"Unknown\n10")

if __name__ == "__main__":
    # Check if the filename is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    data = read_sample_data(filename)

    # Classify the population based on the provided conditions
    classify_population(data)
