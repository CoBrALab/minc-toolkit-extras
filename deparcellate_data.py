import argparse
import numpy as np
import os
import sys



def deparcellate_data(data, parcellation, output, labels=None):
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

    np.savetxt(output, vertexwise_data)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Maps region-averaged values to their corresponding surface vertices'
    )

    parser.add_argument('input', help='path to .txt file of region-averaged values', type=str)
    parser.add_argument('parcellation', help='path to .txt file of labels for each vertex', type=str)
    parser.add_argument('output', help='path to .txt file to store vertex-wise map of regional values', type=str)
    parser.add_argument('--labels', help='optionally specifiy labels for each region in the input file, otherwise assumes regions are in ascending numeric order (according to values in the parcellation file)', default=None, type=str)

    args = parser.parse_args()

    deparcellate_data(data=args.input, parcellation=args.parcellation, output=args.output, labels=args.labels)