"""
Meteorological functions for wind profile
interpolation

V0 - 26.02.2017

TODO: 1.) Unstable Strafitication for Logprof Correction
      2.) Interpolation for Ekman Layer

Dependencies: numpy
"""
import numpy as np

KAPPA=0.4

def logprof(z,ustar,znot,L=np.inf,d=0):
    # Logarithmic Wind Profiles
    # IN: z (des. height) [m], ustar (friction velocity) [m/s],
    #     znot (roughness length) [m], d (displacement height) [m],
    #     L (Obhukov Length) [m]
    a=0.5
    A=1.0
    B=2.0/3.0
    C=5.0
    D=0.35
    # Calculate Correction Functions
    if L==np.inf or L==-1.0*np.inf: # Neutral Case
        psi=0
    else:
        if L>0: # Stable Case
            if ((z/L)>0 and (z/L)<=0.5):
                psi=-1.0*a*z/L
            elif ((z/L)>0.5 and (z/L)<=7):
                psi=(A*z/L)+(B*(z/L-C/D))
            else:
                psi=np.exp(-D*z/L)+(B*C/D)
        elif L<=0: # Unstable Case
            psi=0
    # Calculate Wind Speed Profile with Correction
    uint=(ustar/KAPPA)*(np.ln(((z-d)/znot))-psi)
    return(uint)


def alphaprof(z,uref,zref,alpha):
    # Empirical Power Law Profile
    # IN: z (des. Height) [m], Uref [m/s] at zref [m], shear exponent (alpha)
    uint=uref*((z/zref)**(alpha))
    return(uint)
