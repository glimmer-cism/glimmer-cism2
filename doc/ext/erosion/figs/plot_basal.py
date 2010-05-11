#! /usr/bin/env python
# plot basal sediment models

import Numeric, PyGMT, math
from  pygsl import integrate

def calc_alpha(Nz,phi):
    return (1/(Nz*math.tan(math.radians(phi))))
def calc_beta(N0,Nz,c,phi):
    return - (N0+c/math.tan(math.radians(phi)))/Nz

def calc_sigma(z,N0,Nz,c,phi):
    """Calculate yield strength.

    z: depth
    N0: effective pressure
    Nz: effective pressure gradient
    c: cohesion
    phi: angle of internal friction"""

    return (N0+Nz*z)*math.tan(math.radians(phi))+c

def calc_za(tau,N0,Nz,c,phi):
    """Calculate depth of deforming layer.

    tau: basal shear stress
    N0: effective pressure
    Nz: effective pressure gradient
    c: cohesion
    phi: angle of internal friction"""

    za = Numeric.zeros(len(tau),Numeric.Float)
    beta  = calc_beta(N0,Nz,c,phi)
    alpha = calc_alpha(Nz,phi)
    
    za = Numeric.minimum(0., alpha*tau+beta)
    
    return za

def flow1(z,params):
    """flow law 1.

    tau: basal shear stress
    z: depth
    A: flow law factor
    m: exponent of effective pressure
    n: exponent of basal shear stress
    N0: effective pressure
    Nz: effective pressure gradient"""

    tau = params[0]
    A   = params[1]
    m   = params[2]
    n   = params[3]
    N0  = params[4]
    Nz  = params[5]
    c   = params[6]
    phi = params[7]
    upper  = params[8]

    if za == None:
        return A*tau**n/(N0+Nz*z)**m
    else:
        return (upper-z)*A*tau**n/(N0+Nz*z)**m

def flow2(z,params):
    """flow law 1.

    tau: basal shear stress
    z: depth
    A: flow law factor
    m: exponent of effective pressure
    n: exponent of basal shear stress
    N0: effective pressure
    Nz: effective pressure gradient
    c: cohesion
    phi: angle of internal friction"""

    tau = params[0]
    A   = params[1]
    m   = params[2]
    n   = params[3]
    N0  = params[4]
    Nz  = params[5]
    c   = params[6]
    phi = params[7]
    upper  = params[8]

    if upper == None:
        return A*abs(tau-calc_sigma(z,N0,Nz,c,phi))**n/(N0+Nz*z)**m
    else:
        return (upper-z)*A*abs(tau-calc_sigma(z,N0,Nz,c,phi))**n/(N0+Nz*z)**m
     
def calc_vsld(tau,A,m,n,N0,Nz,c,phi,sigma=False):

    za = calc_za(tau,N0,Nz,c,phi)
    v = Numeric.zeros(len(tau),Numeric.Float)

    if sigma:
        f = flow2
    else:
        f = flow1

    for i in range(0,len(v)):
        sys = integrate.gsl_function(f,[tau[i],A,m,n,N0,Nz,c,phi,None])
        flag,result,error,num = integrate.qng(sys,za[i],0.,1e-8, 1e-8)
        if flag != 0:
            print flag
        v[i] = result
    return v

def calc_vavg(tau,A,m,n,N0,Nz,c,phi,sigma=False):
    """Calculating average sediment velocities.
    
    tau: basal shear stress
    A: flow law factor
    m: exponent of effective pressure
    n: exponent of basal shear stress
    N0: effective pressure
    Nz: effective pressure gradient
    c: cohesion
    phi: angle of internal friction"""

    za = calc_za(tau,N0,Nz,c,phi)
    v = Numeric.zeros(len(tau),Numeric.Float)
    
    if sigma:
        f = flow2
    else:
        f = flow1

    for i in range(0,len(v)):
        sys = integrate.gsl_function(f,[tau[i],A,m,n,N0,Nz,c,phi,0.])
        flag,result,error,num = integrate.qng(sys,za[i],0.,1e-6, 1e-6)
        if flag != 0:
            print flag
        if za[i]!=0:
            v[i] = -result/za[i]
    return v    

if __name__ == '__main__':

    c   = [3.75,15,30,70,70,150,150,250]
    phi = [32,30,27,32,30,32,32,35]


    plot = PyGMT.Canvas('plot_basal.ps',size='A4')
    plot.defaults['LABEL_FONT_SIZE']='12p'
    plot.defaults['ANOT_FONT_SIZE']='10p'

    ay = 4.
    dy = 0.5

    leg = 4

    bigarea = PyGMT.AreaXY(plot,size=[30,30],pos=[0,0])


    tau = Numeric.arange(0.,100.,1.,Numeric.Float)

    key = PyGMT.KeyArea(bigarea,size=[15,2.5])
    key.num = [2,4]

    area_za = PyGMT.AutoXY(bigarea,size=[15,ay],pos=[0,leg])
    area_za.xlabel="shear stress [kPa]"
    area_za.ylabel="depth of deforming sediment layer [m]"
    area_za.axis = 'WeSn'

    area_vsld = PyGMT.AutoXY(bigarea,size=[15,ay],pos=[0,leg+ay+dy])
    area_vsld.ylabel="sliding velocity [m/a]"
    area_vsld.axis = 'Wesn'

    area_vavg = PyGMT.AutoXY(bigarea,size=[15,ay],pos=[0,leg+2*(ay+dy)])
    area_vavg.ylabel="average sediment velocity [m/a]"
    area_vavg.axis = 'Wesn'

    def p(tau,N0,Nz,c,phi, colour, name):
        key.plot_line(name,'-W2/%s'%colour)
        
        area_za.line('-W2/%s'%colour,tau,calc_za(tau, N0,Nz,c,phi))
        #area_vsld.line('-W2/%s'%colour,tau,calc_vsld(tau, 34.8,1.8,1.33, N0,Nz,c,phi,False))
        area_vsld.line('-W2/%st'%colour,tau,calc_vsld(tau,107.11,1.35,0.77 ,N0,Nz,c,phi,True))
        area_vsld.line('-W2/%sta'%colour,tau,calc_vsld(tau, 380.86,2,1,N0,Nz,c,phi,True))

        #area_vavg.line('-W2/%s'%colour,tau,calc_vavg(tau, 34.8,1.8,1.33, N0,Nz,c,phi,False))
        area_vavg.line('-W2/%s'%colour,tau,calc_vavg(tau,107.11,1.35,0.77 ,N0,Nz,c,phi,True))
        area_vavg.line('-W2/%sta'%colour,tau,calc_vavg(tau, 380.86,2,1,N0,Nz,c,phi,True))


    # Breidamerkurkoekull
    p(tau,50., -10., 3.75, 32., '255/0/0', "Breidamerkurkoekull, N=50kPa")

    # typical till
    p(tau,50., -10., 15, 30, '0/255/0', "typical till, N=50kPa")

    # Breidamerkurkoekull
    p(tau,20., -10., 3.75, 32., '0/0/255', "Breidamerkurkoekull, N=20kPa")

    # typical till
    p(tau,20., -10., 15, 30, '0/255/255', "typical till, N=20kPa")
    
    area_za.coordsystem()
    area_vsld.coordsystem()
    area_vavg.coordsystem()
    
    plot.close()

