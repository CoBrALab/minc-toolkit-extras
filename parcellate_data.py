import argparse
import numpy as np
import os
import sys



def parcellate_data(data, parcellation, output, mask=None, save_labels=True):
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
    # Then average values in each parcel, store as a list, and save
    labels = [label for label in np.unique(parcellation) if label > 0]
    parcellated_data = [np.average(data[np.where(parcellation==label)]) for label in labels]
    np.savetxt(output, parcellated_data)


    # LABELS: if we want to write the labels as well, save them in the same directory as the out_file
    if save_labels:
        # If out_file is at the end of a provided path, write it there
        out_dir = output.rsplit('/',1)[0]

        # Otherwise, write it to current directory
        if not os.path.isdir(out_dir):
            out_dir = os.getcwd()
        
        np.savetxt(f'{out_dir}/parcel_labels.txt', labels)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Parcellates surface data and outputs a .txt file with region-averaged values'
    )

    parser.add_argument('input', help='path to .txt file of vertex-wise values to be parcellated', type=str)
    parser.add_argument('parcellation', help='path to .txt file of labels for each vertex', type=str)
    parser.add_argument('output', help='path to .txt file to store region-averaged values', type=str)
    parser.add_argument('--mask', help='path to .txt file indicating included vertices with a 1 and excluded vertices with a 0', type=str)
    parser.add_argument('--save_labels', help='specify if labels for each region should be saved, on by default', dest='save_labels', default=True, action='store_true')
    parser.add_argument('--no_save_labels', help='turn off label saving feature, if desired', dest='save_labels', action='store_false')

    args = parser.parse_args()

    parcellate_data(data=args.input, parcellation=args.parcellation, output=args.output, mask=args.mask, save_labels=args.save_labels)