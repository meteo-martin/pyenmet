"""
Very basic meteorological conversion functions

V0 - 26.02.2017
"""

RLAIR=287.058
CPAIR=1005
KAPPA=RLAIR/CPAIR

def rhoair(tink,pinpa):
    # Calculates the density of dry air
    # IN: Temperature [Kelvin], Pressure [Pa]
    densair=pinpa/(RLAIR*tink)
    return(densair)

def pottemp(tink,pinpa):
    # Calculate potential temperature
    # IN: Pressure [Pa], Temperature [Kelvin]
    pref=1000.0*100.0
    pot=tink*((pref/pinpa)**(KAPPA))
    return(pot)

def magnus(tink):
    # Magnus formula for calculating saturation vapour pressure (in Pa)
    # from Etling, 2008
    # IN: Temperature [Kelvin]
    esp=6.10*10.0**((7.45*(tink-273.0))/(tink-38.0))
    esp=esp*100.0
    return(esp)

def qfromrh(rh,tink,pinpa):
    # Calculates the specific humidity
    # IN: Relative Humidity [Percent], Temperature [Kelvin], Pressure [Pa]
    vapp=(rh/100.0)*(magnus(tink))
    q=0.622*(vapp/pinpa)
    return(q)

def vtemp(rh,tink,pinpa):
    # Calculates the virtual temperature
    # IN: Relative Humidity [Percent],  Temperature [Kelvin], Pressure [Pa]
    vtemp=(1.0+0.608*qfromrh(rh,tink,pinpa))*tink
    return(vtemp)

def vpottemp(ptemp,rh,tink,pinpa):
    # Calculates Virtual Potential Temperature
    # IN: Potential Temperature [Kelvin], Relative Humidity [Percent]
    #     Temperature [Kelvin], Pressure [Pa]
    vpot=(1.0+0.608*qfromrh(rh,tink,pinpa))*ptemp
    return(vpot)
