"""
Meteorological functions for assessment of atmospheric
stratification

V0 - 26.02.2017

Dependencies: numpy
"""
import numpy as np

GRAVITY=9.81

def ribulk(pot2,pz2,pot1,pz1,u2,uz2,u1,uz1):
    # Calculate Richardson Bulk Number
    # IN: Potten
    potm=(pot2/pot1)/2.0
    dudz=(u2-u1)/(uz2-uz1)
    dTdz=(pot2-pot1)/(pz2-pz1)
    rib=(potm/GRAVITY)*(dTdz/(dudz**2.0))
    return(rib)

def shearexp(u2,z2,u1,z1):
    # Shear/Power Law Exponent (alpha)
    alpha=(np.log(u2)-np.log(u1))/(np.log(z2)-np.log(z1))
    return(alpha)
