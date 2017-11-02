"""
Very basic conversion functions for
meteorological quantities

V0 - 26.02.2017 - Thermodynamic Variables
V1 - 02.11.2017 - added wind vector

Dependencies: numpy
"""
import numpy as np

# Meteorological Constants
RLAIR=287.058
CPAIR=1005.0
GRAV=9.81
KAPPA=RLAIR/CPAIR


def uv2ws(u,v):
    """ Calculates the horizontal wind speed
    from the wind vector
    -------------------------------------------------
    In:  U - Longitudinal Wind Vector Component [m/s]
         V - Lateral Wind Vector Component [m/s]
    -------------------------------------------------
    Out: Horizontal Wind Speed [m/s]
    """
    ws=np.sqrt(u*u+v*v)
    return(ws)


def uv2wd(u,v):
    """ Calculates the horizontal wind directon
    from the wind vector
    -------------------------------------------------
    In:  u - Longitudinal Wind Vector Component [m/s]
         v - Lateral Wind Vector Component [m/s]
    -------------------------------------------------
    Out: wd - wind direction [degrees]
    """
    wsabs=uv2ws(u,v)
    wdrad=np.degrees(np.arctan2(u/wsabs,v/wsabs))
    wd=np.mod((wdrad+180.0),360.0)
    return(wd)


def degC2K(tc):
    """ Converts degrees Celsius to Kelvin
    ----------------------------------------
    In:  tc - Temperature [degrees C]
    ----------------------------------------
    Out: tk - Temperature [Kelvin]
    """
    tk=tc+273.15
    return(tk)


def hPa2Pa(phpa):
    """ Converts pressure in hPa to Pacal
    ----------------------------------------
    In:  tc - Temperature [degrees C]
    ----------------------------------------
    Out: tk - Temperature [Kelvin]
    """
    pinpa=phpa*100.0
    return(pinpa)


def rhoair(tc,p):
    """ Calculates the density of dry air
    ------------------------------------------
    In: tc - Temperature [degree C]
         p - Pressure [hPa]
    ------------------------------------------
    Out: dens - air density [kg/m^3]
    """
    tink=degC2K(tc)
    pinpa=hPa2Pa(p)
    densair=pinpa/(RLAIR*tink)
    return(densair)


def pottemp(tc,p):
    """ Calculates potential temperature
    ----------------------------------------
    In: tc - Temperature [degree C]
         p - Pressure [hPa]
    ----------------------------------------
    Out: pot - potential temperature [K]
    """
    PREF=1000.0*100.0
    tink=degC2K(tc)
    pinpa=hPa2Pa(p)
    pot=tink*((PREF/pinpa)**(KAPPA))
    return(pot)


def magnus(tc):
    """ Magnus formula for calculating saturation
        vapour pressure (in Pa) from Etling, 2008
    -----------------------------------------------------
    In:  tink - Temperature [Kelvin]
    -----------------------------------------------------
    Out: esp - Saturation vapour pressure [Pa]
    """
    tink=degC2K(tc)
    esp=6.10*10.0**((7.45*(tink-273.0))/(tink-38.0))
    esp=esp*100.0
    return(esp)


def rh2q(rh,tc,p):
    """ Calculates specific humidity from relative humidity
    ----------------------------------------
    In: rh - Relative humidity [percent]
        tc - Temperature [degree C]
         p - Pressure [hPa]
    ----------------------------------------
    Out: q - specific humidity [kg/kg]
    """
    tink=degC2K(tc)
    pinpa=hPa2Pa(p)
    vapp=(rh/100.0)*(magnus(tink))
    q=0.622*(vapp/pinpa)
    return(q)


def vtemp(rh,tc,p):
    """ Calculates the virtual temperature
    -----------------------------------------------------
    In: rh - relative humidity [Percent]
        tc - temperature [degree C]
         p - pressure [hPa]
    -----------------------------------------------------
    Out: vtemp - virtual temperature [K]
    """
    tink=degC2K(tc)
    vtemp=(1.0+0.608*rh2q(rh,tc,p))*tink
    return(vtemp)


def vpottemp(rh,tc,p):
    """ Calculates the virtual potential temperature
    -----------------------------------------------------
    In: rh - relative humidity [Percent]
        tc - temperature [degree C]
         p - pressure [hPa]
    -----------------------------------------------------
    Out: vpot - virtual potential temperature [K]
    """
    vpot=(1.0+0.608*rh2q(rh,tc,p))*pottemp(tc,p)
    return(vpot)


