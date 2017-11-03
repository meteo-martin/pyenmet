"""
Collection of Wake Models

TODO: Check Larsen Model

"""

import numpy as np
import matplotlib.pyplot as plt

def jensen(u0,D,ct,k,x):
    """ Jensen wake model
    --------------------------------
    In: u0 - inflow wind speed [m/s]
         D - rotor diameter [m]
        ct - thrust coefficient [-]
         k - wake coefficient [-]
    --------------------------------
    Out: vx - wind speed at Distance x [m]
         dw - wake width [m]
    """
    if ct>1.0:
        ct=1.0
    dwf=2.0*k*x # Linear Expanding Wake Width factor
    dw=dwf+D
    ww=(dwf/D)
    vx=u0*(1.0-((1.0-np.sqrt(1.0-ct))/(1+ww**2.0)))
    return([vx,dw]) # Returns Width of Wake and Velocity

def larsen(u0,D,H,ct,Ia,r,x,x0):
    """ Jensen wake model
    according to Master Thesis by J. Renkema
    --------------------------------
    In: u0 - inflow wind speed [m/s]
         D - rotor diameter [m]
        ct - thrust coefficient [-]
        Ia - ambient turbulence intensity [percent]
         r - radial distance [m]
         x - lateral distance [x]
    --------------------------------
    Out: vx - wind speed at position x,r [m]
    """
    Ia=Ia/100.0

    # Rotor Area
    a=np.pi*((D/2.0)**2.0)
    # Effective Rotor Diameter
    deff=D*np.sqrt((1.0+np.sqrt(1.0-ct))/(2.0*np.sqrt(1.0-ct)))

    rnb=np.max([1.08*D,(1.08*D+21.7*D*(Ia-0.05))])

    # Wake Radius at a distance of 9.5 D
    r95=0.5*(rnb+np.min([H,rnb]))

    # Rotor Psition with respect to applied Coordinated System
    x0=(9.5*D)/(((2.0*r95/deff)**3.0)-1.0) # position of the rotor

    # Constant Related to the Prandtl Mixing Length
    c0=((deff/2.0)**(5.0/2.0))
    c1=c0*((105.0/(2.0*np.pi))**(-1.0/2.0))*((ct*a*x0)**(-5.0/6.0))

    # Calculate Wind Speed Deficit
    f0=(ct*a*(x+x0)**(-2.0))**(1.0/3.0) #OK
    f1=-1.0*(u0/9.0)*f0 #OK
    f2=(3.0*(c1**2.0)*ct*a*(x+x0))**(-1.0/2.0) #OK
    f3=(35.0/2.0*np.pi)**(3.0/10.0) #OK
    f4=(3.0*c1**2.0)**(-1.0/5.0) #OK
    vwx=-1.0*((f1*(r**(3.0/2.0)*f2-f3*f4)**2.0)-u0)

    # Calculate Wake Radius
    g1=(35.0/(2.0*np.pi))**(1.0/5.0)
    g2=(3.0*c1**2.0)**(1.0/5.0)
    g3=(ct*a*(x+x0))**(1.0/3.0)
    rw=g1*g2*g3

    return([vx,rw])


def main():
    pass

    # Input Parameters to the Model:
    ct=0.8
    uin=8.0
    D=126.0
    H=130.0
    k=0.05

    dist=np.arange(0,10000.0)
    [vj,wwj]=jensen(uin,D,ct,k,dist)
    [vl,wwl]=larsen(uin,D,H,ct,10.0,0.0,dist,0.0)

    # Plots
    f,ax=plt.subplots(2,sharex=True)
    ax[0].plot(dist/D,vj/uin,label="Jensen")
    ax[0].plot(dist/D,vl,label="Larsen")

    ax[0].set_ylabel("Wind Speed $U$/$U_{in}$ [-]")
    ax[1].plot(dist/D,wwj/D,label="Jensen")
    ax[1].plot(dist/D,wwl/D,label="Larsen")

    ax[1].set_ylabel("Wake Width $Wake_{width}/D$ [-]")
    plt.xlabel("Distance [$D$]")
    plt.legend(loc=4)
    plt.show()

if __name__ == "__main__":
    main()
