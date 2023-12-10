#!/usr/bin/env python

import argparse
import numpy as np
import os
import sys



def parcellate_data(data, parcellation, mask=None):
    # INPUT: Load in data and parcellation files (assuming .txt for now)
    data = np.genfromtxt(data)
    parcellation = np.genfromtxt(parcellation)

    if len(data) != len(parcellation):
        sys.exit('Error: Lengths of input data and parcellation must be identical')


    # MASKING: if mask file provided, remove any excluded vertices from data and parcellation
    if mask is not None:
        mask = np.genfromtxt(mask)
        
        if len(data) != len(mask):
            sys.exit('Error: Lengths of input data/parcellation and mask must be identical')

        if not np.all((mask==0) | (mask==1)):
            sys.exit('Error: Mask must be binary (0/1)')

        data = data[mask==1]
        parcellation = parcellation[mask==1]


    # PARCELLATION: get unique parcel labels - ignoring anything negative or equal to 0
    # Then average values in each parcel, store as a list, and return
    labels = [label for label in np.unique(parcellation) if label > 0]
    parcellated_data = [np.average(data[np.where(parcellation==label)]) for label in labels]

    return parcellated_data, labels



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generic function for parcellating CIVET surface measures'
    )

    parser.add_argument('input', help='vertex-wise input data to be parcellated. Format: path to .txt file with one value for each vertex', type=str)
    parser.add_argument('parcellation', help='parcellation scheme. Format: path to .txt file with one label for each vertex', type=str)
    parser.add_argument('output', help='output file to store region-averaged values. Format: path to .txt file', type=str)
    parser.add_argument('--mask', help='optional mask to exclude specified vertices from the regional averaging. Format: path to .txt file marking included vertices with a 1 and excluded vertices with a 0', type=str)
    parser.add_argument('--label-output', help='optional output file containing labels for each region-averaged value in output. Format: path to .txt file', type=str)

    args = parser.parse_args()
    
    # Get parcellated data and labels for each parcel
    parcellated_data, labels = parcellate_data(data=args.input, parcellation=args.parcellation, mask=args.mask)
    
    # Save parcellated data, optionally save labels as well
    np.savetxt(args.output, parcellated_data)
    
    if args.label_output is not None:
        np.savetxt(args.label_output, labels)
        