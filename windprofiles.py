"""
Meteorological functions for wind profile
interpolation

V0 - 26.02.2017

TODO: 1.) Interpolation for Ekman Layer

Dependencies: numpy
"""
import numpy as np

KAPPA=0.4

def logprof(z,ustar,znot,L=np.inf, disph=0.0):
    """
    Logarithmic Wind Profile
    IN: z (des. height) [m], ustar (friction velocity) [m/s],
        znot (roughness length) [m], d (displacement height) [m],
        L (Obhukov Length) [m]
    """
    # Calculate Correction Function psi (from Emeis 2013)
    if L==np.inf or L==-1.0*np.inf: # Neutral Case
        psi=0.0
    else:
        if L>0.0: # Stable Case
            # Constants:
            asm=0.5
            A=1.0
            B=2.0/3.0
            C=5.0
            D=0.35
            # Correction
            if ((z/L)>0.0 and (z/L)<=0.5):
                psi=-1.0*asm*z/L
            elif ((z/L)>0.5 and (z/L)<=7.0):
                psi=(A*z/L)+(B*(z/L-C/D))
            else:
                psi=np.exp(-D*z/L)+(B*C/D)
        elif L<=0.0: # Unstable Case
            # Constants:
            b=16.0
            x=(1.0-(b*z/L))**0.25
            x1=(2.0*np.log((1.0+x)/2.0))
            x2=np.log((1.0+x**2.0)/2.0)
            x3=2.0*np.arctan(x)
            psi=x1+x2-x3+(np.pi/2.0)
    # Calculate Wind Speed Profile with Correction
    uint=(ustar/KAPPA)*(np.log((z-disph)/znot)-psi)
    return(uint)

def alphaprof(z,uref,zref,alpha):
    """
    Empirical Power Law Profile
    IN: z (des. Height) [m], Uref [m/s] at zref [m], shear exponent (alpha)
    """
    uint=uref*((z/zref)**(alpha))
    return(uint)
