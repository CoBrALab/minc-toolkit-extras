#!/usr/bin/env python

import argparse
import numpy as np
import os
import sys



def deparcellate_data(data, parcellation, labels=None):
    # INPUT: Load in data and parcellation files (assuming .txt for now)
    data = np.genfromtxt(data)
    parcellation = np.genfromtxt(parcellation)

    # Either load in labels if specified, or get them from the parcellation
    # Makes sure all labels are positive
    if labels is None:
        labels = [label for label in np.unique(parcellation) if label > 0]
    else:
        labels = np.genfromtxt(labels)

        if not np.all(labels > 0):
            sys.exit('Error: All parcel labels must be positive')

    if len(data) != len(labels):
        sys.exit('Error: Must be one value for each parcel label')



    # DEPARCELLATION: pre-define vertex-wise array
    # For each label, get corresponding value, fill vertexwise array accordingly
    # Then save output
    vertexwise_data = np.zeros(np.shape(parcellation))

    for (idx, label) in enumerate(labels):
        parcel_value = data[idx]
        to_fill = np.where(parcellation == label)
        vertexwise_data[to_fill] = parcel_value

    return vertexwise_data




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generic function for mapping parcellated CIVET surface measures to their corresponding vertices'
    )

    parser.add_argument('input', help='region-averaged data. Format: path to .txt file with one value for each region', type=str)
    parser.add_argument('parcellation', help='parcellation scheme. Format: path to .txt file with one label for each vertex', type=str)
    parser.add_argument('output', help='output file to store region-averaged values at their corresponding vertices. Format: path to .txt file', type=str)
    parser.add_argument('--labels', help='optionally specify labels for each region in the input data. If not specified, assumes regions are in ascending numeric order, based on the unique values in parcellation. Format: path to .txt file with one label for each region', default=None, type=str)

    args = parser.parse_args()

    # Get vertex-wise data, save
    vertexwise_data = deparcellate_data(data=args.input, parcellation=args.parcellation, labels=args.labels)
    np.savetxt(args.output, vertexwise_data)
    