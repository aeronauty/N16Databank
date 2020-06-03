#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 14:46:35 2020

@author: harrysmith

This module imports and interprets the data contained within N16TABLE.dat and
creates spline interpolant functions for the user inputs - there are two useful
functions:
    MakeInterpolants(Cld_in, tc_in): returns 4 interpolant functions
    - Cl_interpolant_subsonic, Cl_interpolant_supersonic, Cd_interpolant_subsonic, Cd_interpolant_supersonic
    
    ClCd(Cld_in, tc_in, mach_in, alpha_in): runs MakeInterpolants and then
    determines the Cl and Cd for the used inputs given. This is good for testing, but not recommended to use
    as part of any code
    
The data were taken from wind tunnel tests of untracable legacy but have shown
good correlation to propeller performance data.

The alpha range is reduced for the supersonic range hence the separation of interpolants
"""

import numpy as np
import csv

# Load the databank in from the .dat file
datafile = open('N16TABLE.dat', 'r')
datareader = csv.reader(datafile)
databanks = [] # Empty list to put our databanks in
for row in datareader:
    # print(row)
    temprow = row[0].split()
    if len(temprow) == 1: # This condition is only met if at a header row
        tablnum = int(temprow[0]) # Get the value of the databank (the single nu)
        irow = 0
        
        if tablnum > 1:
            databanks.append(databank)

        databank = np.empty((20,9)) # Initialise an empty databank to store the next one in             
        irow = 0
    else:
        tablerow = np.array(temprow)
        databank[irow, :] = tablerow[:-1]
        irow = irow + 1
        
databanks.append(databank) # Puts the last entry in the databank

# What we have at this point is 312 sets of tables - each of which contains the Cl and Cd for a given alpha and mach
# with the first ten rows being Cl and the second ten rows being Cd
#
# The ten rows of each corresponds to the t_c ratio, and the columns refer to the design lift coefficients.
# The tested values of these are given below

cldes = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
tmonc = np.array([0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3])

""" 
For a given Cld and and t/c ratio, each of the 312 x 2 tables should be interpolated

This will yield the full alpha/mach variation of Cl and Cd for a given aerofoil
"""      

def MakeInterpolants(Cld_in=.2, tc_in=.2):
    """Creates 4 interpolant objects for the design conditions given

    Sub/supersonic Cl and Cd interpolants - f(alpha [deg], Mach)
    """
    # Logic to handle the replacement of values outside of range
    if Cld_in < cldes.min():
        print("Design lift coefficient in ({:1.2f}) outside of table range - setting to Cld = {:1.2f}".format(Cld_in, cldes.min()))
        Cld_in = cldes.min()
    elif Cld_in > cldes.max():
        print("Design lift coefficient in ({:1.2f}) outside of table range - setting to Cld = {:1.2f}".format(Cld_in, cldes.max()))
        Cld_in = cldes.max()   
        
    if tc_in < tmonc.min():
        print("Thickness ratio in ({:1.2f}) outside of table range - setting to t/c = {:1.2f}".format(tc_in, tmonc.min()))
        tc_in = tmonc.min()
    elif Cld_in > cldes.max():
        print("Thickness ratio in ({:1.2f}) outside of table range - setting to t/c = {:1.2f}".format(tc_in, tmonc.min()))
        tc_in = tmonc.max()   
    
    # Make some empty arrays to put the 312 values of lift and drag into
    Cls = np.zeros((1, 312))
    Cds = np.zeros((1, 312))
    
    # Do the interpolation
    from scipy.interpolate import interp2d
    for i in range(0, 312):
        Cltable = databanks[i][:10, :]
        Clinterpolant = interp2d(cldes, tmonc, Cltable, kind='linear')
        Cls[0, i] = Clinterpolant(Cld_in, tc_in)
        
        Cdtable = databanks[i][10:, :]
        Cdinterpolant = interp2d(cldes, tmonc, Cdtable, kind='linear')
        Cds[0, i] = Cdinterpolant(Cld_in, tc_in)
    
    # These output data need to be separated into the two alpha ranges - we'll 
    # call them subsonic/supersonic
    
    # Alpha and Mach values in:
    alpha =[-6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36]        
    mach = [.3,.4123,.5,.5745,.6403,.7,.755,.8062,.8544,.9,.9434,.9849,1.025,1.063,1.1,1.136,1.368,1.6]
    
    ind_boundary = 22*12 
    
    Cl_subsonic = Cls[0, :ind_boundary]
    Cd_subsonic = Cds[0, :ind_boundary]
    
    Cl_supersonic = Cls[0, ind_boundary:]
    Cd_supersonic = Cds[0, ind_boundary:]
    
    # Reshape this into an array with M as the row and alpha as the columns
    Cl_subsonic = Cl_subsonic.reshape((12, 22))
    Cd_subsonic = Cd_subsonic.reshape((12, 22))
    Cl_supersonic = Cl_supersonic.reshape((6, 8))
    Cd_supersonic = Cd_supersonic.reshape((6, 8))
    
    # Make some interpolants for supersonic and subsonic
    Cl_interpolant_subsonic = interp2d(alpha, mach[:12], Cl_subsonic, kind='cubic')
    Cl_interpolant_supersonic = interp2d(alpha[:8], mach[12:], Cl_supersonic, kind='cubic')
    Cd_interpolant_subsonic = interp2d(alpha, mach[:12], Cd_subsonic, kind='cubic')
    Cd_interpolant_supersonic = interp2d(alpha[:8], mach[12:], Cd_supersonic, kind='cubic')
    
    return Cl_interpolant_subsonic, Cl_interpolant_supersonic, Cd_interpolant_subsonic, Cd_interpolant_supersonic

def ClCd(Cld_in=.2, tc_in=.2, mach_in=0.3, alpha_in=4):
    """Creates the interpolants for Cl and Cd and then determines the values based on input
    
    This method is less efficient to use in a program as it has to create every interpolant for each table each time the program is run. Good for a single case/debugging.
    """
    Cl_interpolant_subsonic, Cl_interpolant_supersonic, Cd_interpolant_subsonic, Cd_interpolant_supersonic = MakeInterpolants(Cld_in, tc_in)

    # Do the interpolation
    if mach_in <= 1:
        Cl = Cl_interpolant_subsonic(alpha_in, mach_in)
        Cd = Cd_interpolant_subsonic(alpha_in, mach_in)
    else:
        Cl = Cl_interpolant_supersonic(alpha_in, mach_in)
        Cd = Cd_interpolant_supersonic(alpha_in, mach_in)
    return Cl, Cd
