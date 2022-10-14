#!/usr/bin/env python

"""@package docstring
File: objects.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import numpy as np
import math
import yaml
import vtk
from numba import jit, vectorize
import argparse


class filament():

    def __init__(self, line):
        if type(line) == bytes:
            line = line.decode()
        self.info = line.split()
        self.fil_type = str(self.info[0])
        self.gid = int(self.info[1])

        dat = np.asarray(self.info[2:], dtype=np.double)
        self.radius = dat[0]
        self.minus_end = dat[1:4]
        self.plus_end = dat[4:7]
        self.vec = self.plus_end - self.minus_end

    def parse(self):
        dat = np.asarray(self.info[2:8], dtype=np.double)
        self.lengthxy = np.linalg.norm(self.vec[:-1])
        self.length = np.linalg.norm(self.vec)
        if self.length != 0:
            self.orientation = self.vec / self.length
            self.theta = -np.arctan2(self.orientation[1], self.orientation[0])
        else:
            self.orientation = self.vec
            self.theta = 0

    def get_com(self):
        return .5 * (self.plus_end + self.minus_end)

    def get_dat(self):
        return self.info[1:]


class protein():

    def __init__(self, line):
        if type(line) == bytes:
            line = line.decode()
        self.info = line.split()
        self.xlp_type = str(self.info[0])
        self.gid = int(self.info[1])

    # def parse(self):
    #     dat = np.asarray(self.info[2:8], dtype=np.double)
    #     self.radius = dat[0]
    #     self.minus_end = dat[1:4]
    #     self.plus_end = dat[4:7]

    #     self.vec = self.plus_end - self.minus_end
    #     self.lengthxy = np.linalg.norm(self.vec[:-1])
    #     self.length = np.linalg.norm(self.vec)
    #     if self.length != 0:
    #         self.orientation = self.vec / self.length
    #         self.theta = -np.arctan2(self.orientation[1], self.orientation[0])
    #     else:
    #         self.orientation = self.vec
    #         self.theta = 0

    # def get_com(self):
    #     return .5 * (self.plus_end + self.minus_end)

    def get_dat(self):
        return self.info[1:]


class con_block():

    """Docstring for con_block. """

    def __init__(self,):
        """ Need only the position of to acting points to start. Other variables
        are added when this is filled."""
        end0 = None
        end1 = None


##########################################
if __name__ == "__main__":
    print("Not implemented.")
