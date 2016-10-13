#!/usr/bin/env python

import kernel_scipy as kernel
import csv
import numpy as np
import automateHeatKernel as ahk
from optparse import OptionParser
import operator
from scipy.sparse import csc_matrix

class EstablishedKernel(kernel.SciPYKernel):

    def __init__(self, kernel_file):
        """
           Input:
                    kernel_file - filename of tab delemited kernel file, as made by kernel_scipy.SciPYKernel
        """
        self.readKernel(kernel_file)


    def readKernel(self,input_file):
        ker=[]
        start=True
        for line in csv.reader(open(input_file,'r'),delimiter='\t'):
            if start:
                self.labels=line[1:]
                start=False
            else:
                ker.append([float(x) for x in line[1:]])

        self.kernel=csc_matrix(np.asarray(ker))
