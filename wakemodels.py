import numpy as np
import matplotlib.pyplot as plt

def jensen(u0,D,ct,k,x):
    # k=Wake Coefficient
    # D=Rotor Diameter
    #
    # Returns: vx - Wind Speed at Distance X
    #
    if ct>1.0:
        ct=1.0
    dwf=2.0*k*x # Linear Expanding Wake Width factor
    dw=dwf+D
    ww=(dwf/D)
    vx=u0*(1.0-((1.0-np.sqrt(1.0-ct))/(1+ww**2.0)))
    return([vx,dw]) # Returns Width of Wake and Velocity

def jensen2field(tpos,u0,vx,ww,xmin,xmax,ymin,ymax,res=10.0,wdir):
    y=np.arange(ymin,ymax,res)
    x=np.arange(xmin,xmax,res)
    vmat=np.empty([len(x),len(y)])
    for i in range(0,len(x)):
        for j in range(0,len(y)):
            if x[i]<tpos[0]:
                vmat[i,j]=u0
            else:
                if np.abs((y[i]-tpos[1]))>ww[]
                    vmat[i,j]=
                else:
                    vmat[i,j]=u0
    return([x,y,vmat])


def larsen(u0,D,ct,c1,r,x):
    a=np.pi*((D/2.0)**2.0) # Rotor Area

    # Related to
    c1=(deff/2.0)**5.0/2.0*(105.0/2.0*np.pi)**(-1.0/2.0)*(ct-a*x0)**(-5.0/6.0)
    # position of rotor
    x0=9.5*D

    f0=(ct*a*x**(-2.0))**(1.0/3.0)
    f1=-1.0*(u0/9.0)*f0
    f2=(3.0*(c1**2.0)*ct*a*x)**(-1.0/2.0)
    f3=(35.0/2.0*np.pi)**(3.0/10.0)
    f4=(3.0*c1**2.0)**(-1.0/5.0)
    vx=-1.0*((f1*(r**(3.0/2.0)*f2-f3*f4)**2.0)-u0)
    return(vx)


def main():
    pass

    # Input Parameters to the Model:
    ct=0.8
    uin=8.0
    D=126.0
    k=0.05

    dist=np.arange(0,10000.0)
    [v,wakew]=jensen(uin,D,ct,k,dist)

    # Plots
    f,ax=plt.subplots(2,sharex=True)
    ax[0].plot(dist/D,v/uin)
    ax[0].set_ylabel("Wind Speed $U$/$U_{in}$ [-]")
    ax[1].plot(dist/D,wakew/D)
    ax[1].set_ylabel("Wake Width $Wake_{width}/D$ [-]")
    plt.xlabel("Distance [$D$]")
    plt.show()

if __name__ == "__main__":
    main()
