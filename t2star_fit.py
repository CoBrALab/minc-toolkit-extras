#!/usr/bin/env python

# all the imports go here
#from pyminc.volumes.factory import *
import pyezminc
import sys
import os
import numpy as np
from scipy.optimize import leastsq
from optparse import OptionParser

#def model(coeffs, t, y):
#   return coeffs[0] * np.exp( -t/coeffs[1] ) - y

def model(coeffs, t, y):
   return -t/coeffs[1] + np.log(coeffs[0]) - np.log(y)

if __name__ == "__main__":

    # some argument handling
    usage = 't2star-fit.py "list of echo numbers" mask.mnc <list of echo files> beta.mnc t2star.mnc'
    description = "Description text"
    
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--clobber", dest="clobber",
                      help="clobber output file",
                      type="string")

    (options, args) = parser.parse_args()

    if len(args) < 5:
        parser.error("Incorrect number of arguments")

    outfilename = []
    outfilename.append(args[-1])
    outfilename.append(args[-2])


    #echos= np.array([0.00284, 0.0062, 0.00956, 0.01292, 0.01628, 0.01964, 0.023, 0.02636, 0.02972, 0.03308, 0.03644, 0.0398], dtype=np.float64)
    echos=np.array([float(s) for s in args[0].split()], dtype=np.float64)

    mask = args[1]

    # create an array of volume handles
    input_files = []

    nfiles = len(args)-2
    for i in range(2,nfiles ):
        input_files.append(args[i])
    

    inp=pyezminc.parallel_input_iterator()
    out=pyezminc.parallel_output_iterator()
    inp.open(input_files, mask)
    out.open(outfilename, mask)
    inp.begin()
    out.begin()

    try:
        while True:
            result=np.zeros(shape=[out.dim()],dtype=np.float64,order='C')
            if inp.value_mask():
                guess=np.array([5000, 0.05], dtype=np.float64)
                intensity = inp.value()
                result, flag = leastsq(model, guess, args=(echos, intensity))
            out.value(result)
            
            # move to the next voxel
            inp.next()
            out.next()
            # display something on the screen....
            #sys.stdout.write(".")
            #sys.stdout.flush()

    except StopIteration:
        pass

    del inp
    del out    
