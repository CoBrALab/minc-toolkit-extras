#!/usr/bin/env python
from optparse import OptionParser
import warnings
import SimpleITK as sitk
import numpy as np
from scipy.optimize import curve_fit
import miniutils
import sys

warnings.filterwarnings("error")


def model(t, s0, t2star):
    return s0 * np.exp(-t / t2star)


def do_fit(voxel):
    try:
        if maskdata[voxel] > 0:
            popt, _ = curve_fit(model, echos, inputdata[:, voxel], p0=[350, 0.05])
            return np.array([popt[0], popt[1]], dtype=np.float64)
        else:
            return np.array([0, 0], dtype=np.float64)
    except:
        return np.array([0, 0], dtype=np.float64)


if __name__ == "__main__":

    # some argument handling
    usage = 't2star-fit.py "list of echo numbers" mask.mnc <list of echo files> s0.mnc t2star.mnc'
    description = "Description text"

    parser = OptionParser(usage=usage, description=description)
    parser.add_option(
        "--clobber", dest="clobber", help="clobber output file", type="string"
    )

    (options, args) = parser.parse_args()

    if len(args) < 5:
        parser.error("Incorrect number of arguments")

    s0_outfilename = args[-2]
    t2star_outfilename = args[-1]

    echos = np.array([float(s) for s in args[0].split()], dtype=np.float64)
    mask = args[1]

    # create an array of volume handles
    input_files = args[2:-2]

    if len(input_files) != len(echos):
      sys.exit("Error: length of echos specified not equal to number of input files")

    maskimage = sitk.ReadImage(mask)
    maskdata = sitk.GetArrayFromImage(maskimage)
    img_dim = maskdata.shape
    maskdata = np.ravel(maskdata)

    inputdata = np.zeros(
        (len(echos), img_dim[0] * img_dim[1] * img_dim[2]), dtype=np.float64
    )

    for i, image in enumerate(input_files):
        inputdata[i,] = np.ravel(sitk.GetArrayFromImage(sitk.ReadImage(image)))

    results = miniutils.parallel_progbar(do_fit, range(inputdata.shape[1]))
    results = np.asarray(results)
    s0_est = results[:, 0]
    t2star_est = results[:, 1]

    s0_est = sitk.GetImageFromArray(np.clip(s0_est.reshape(img_dim), 0, 1000))
    s0_est.CopyInformation(maskimage)
    sitk.WriteImage(s0_est, s0_outfilename)

    t2star_est = sitk.GetImageFromArray(np.clip(t2star_est.reshape(img_dim), 0, 1))
    t2star_est.CopyInformation(maskimage)
    sitk.WriteImage(t2star_est, t2star_outfilename)
