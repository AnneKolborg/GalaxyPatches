import numpy as np
import scipy.integrate as integrate
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt

def phi(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    rho_s=p0
    R_s=p1*3.08e21
    Mb=p2*1.989e33
    hb=p3*3.08e21
    Mds=p4*1.989e33
    Mdg=p5*1.989e33
    hr=p6*3.08e21
    hz=p7*3.08e21
    Tmu=p8
    GG=6.67e-8

    rcirc=r*3.08e21
    zcirc=z*3.08e21
    Rspher=np.sqrt(r**2+z**2+1e-10**2)*3.08e21
    
    phi_nfw=-4.0*np.pi*GG*rho_s*R_s**2*np.log(1+Rspher/R_s)/(Rspher/R_s)
    phi_bulge=-GG*Mb/np.sqrt(Rspher**2+hb**2)
    #phi_disc=-GG*(Mds+Mdg)/np.sqrt((rcirc+hr)**2+(np.abs(zcirc)+hz)**2)
    phi_disc=-GG*(Mds+Mdg)/np.sqrt(rcirc**2+(hr+np.sqrt(zcirc**2+hz**2))**2)
    phi_tot = phi_nfw + phi_bulge + phi_disc

    return phi_bulge

def dphi_dr(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    hr=p6*3.08e21
    dh = p6/1.0e4
    out = (phi(r+dh,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)-phi(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8))/(dh*3.08e21)

    return out
    
def dphi_dz(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    hz=p7*3.08e21
    dh = p7/1e4
    out = (phi(r,z+dh,p0,p1,p2,p3,p4,p5,p6,p7,p8)-phi(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8))/(dh*3.08e21)
    
    return out

def vcirc(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    hr=p6*3.08e21
    hz=p7*3.08e21
    Tmu=p8
    dh = p7/1.0e4
    mp=1.66e-24
    kb=1.38e-16
    dpdr=-2*r*1.38e-16*Tmu/1.66e-24/hr**2+dphi_dr(r,dh,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    
    out = (r*3.08e21*dpdr)**0.5

    return out
    
def kappa(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    hr=p6*3.08e21
    dh = p6/1.0e4
    dpdr=dphi_dr(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    dpdr1=dphi_dr(r+dh,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    dpdr0=dphi_dr(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    term1 = 3*np.abs(dpdr)/(r*3.08e21)
    term2 = (np.abs(dpdr1)-np.abs(dpdr0))/(dh*3.08e21)
    out = (term1+term2)**0.5

    return out
    
def rho_scaling(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    rho_s=p0
    R_s=p1*3.08e21
    Mb=p2*1.989e33
    hb=p3*3.08e21
    Mds=p4*1.989e33
    Mdg=p5*1.989e33
    hr=p6*3.08e21
    hz=p7*3.08e21
    Tmu=p8
    GG=6.67e-8
    mp=1.66e-24
    kb=1.38e-16
    dh = p7/1.0e4
    
    rcirc=r*3.08e21
    zcirc=z*3.08e21
    Rspher=np.sqrt(r**2+z**2)*3.08e21

    out = np.exp(-(rcirc/hr)**2)*np.exp(-mp*(phi(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)-phi(r,dh,p0,p1,p2,p3,p4,p5,p6,p7,p8))/kb/Tmu)
        
    return out

def Tmu_norm(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    rho_s=p0
    R_s=p1*3.08e21
    Mb=p2*1.989e33
    hb=p3*3.08e21
    Mds=p4*1.989e33
    Mdg=p5*1.989e33
    hr=p6*3.08e21
    hz=p7*3.08e21
    Tmu=p8
    GG=6.67e-8
    mp=1.66e-24
    kb=1.38e-16
    
    rcirc=r*3.08e21
    zcirc=z*3.08e21
    Rspher=np.sqrt(r**2+z**2+1e-10**2)*3.08e21
    
    #Tmu = hz**2/(zcirc+1.0e-10)*mp/kb/2
    #out = Tmu*dphi_dz(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    out = Tmu
    
    return out

def total_integrand(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    out = 2*np.pi*r*rho_scaling(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    return out

def vertical_integrand(z,r,p0,p1,p2,p3,p4,p5,p6,p7,p8):
    out = rho_scaling(r,z,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    return out

#NFW parameters
h = 0.7 
H0 = h*100*1e5/3.08e24 #in cgs 
M200 = 1e11 # M200 in Msun 
R200 = (2*6.67e-8*1.989e33*M200/200/H0/H0)**(1.0/3.0)/3.08e21 # R200 in kpc
c200 = 8.03*(M200*h/1e12)**(-0.101) # Dutton and Maccio 2014 
R_s = R200/c200 # Scale radius in kpc
nfw_fac = np.log(1+c200)-c200/(1+c200)
rho_s = (M200*1.989e33)/(4.0*np.pi)/(R_s*3.08e21)**3/nfw_fac # Scale density in g/cm3
v200 = np.sqrt(6.67e-8*M200*1.989e33/R200/3.08e21)/1e5
T200 = 3.6e5*(v200/100)**2
print 'NFW Parameters: '
print 'H0 = ',h*100,' km/s/Mpc'
print 'M200 = ',M200,' Msun'
print 'R200 = ',R200,' kpc'
print 'c200 = ',c200
print 'T200 = ',T200
print 'R_s = ',R_s,' kpc'
print 'rho_s = ',rho_s,' g/cm3'
print ' '

#Bulge parameters
Mbulge = 1.0e4  # Bulge mass in Msun
hbulge = 1.0e-1 # Bulge scale radius in kpc
print 'Bulge parameters: '
print 'Mbulge = ',Mbulge,' Msun'
print 'hbulge = ',hbulge,' kpc'
print ' '

#Disc parameters
Mdstar = 5e9 # Stellar disc mass in Msun
Mdgas = 1.0e10 # Stellar disc mass in Msun
hrdisc = 1.6 # Radial scale radius in kpc
hzdisc = 0.8 # Scale height in kpc
print 'Disc parameters: '
print 'Mdstar = ',Mdstar,' Msun'
print 'Mdgas = ',Mdgas,' Msun'
print 'hrdisc = ',hrdisc,' kpc'
print 'hzdisc = ',hzdisc,' kpc'
print ' '

#Compute gas disc normalization and central density
Tmu = 3.0e4
params = (rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
normalization = integrate.nquad(total_integrand, [[0, R200],[-R200, R200]], args=params)
rho0gas = Mdgas*1.989e33/(normalization[0]*3.08e21**3)
print "rho0gas = ", rho0gas,' g/cm3'
print "nH0gas = ", 0.76*rho0gas/1.66e-24,' 1/cm3'
print "T/mu = ",Tmu," K"
print ' '

#Define variables for plots
r = np.linspace(-3,2,60)
z = np.linspace(-3,2,60)
r = 10**r
z = 10**z
Sigmagas=np.zeros((len(r)))

#Compute Sigma_gas profile
i=0
for rr in r:
    params2 = (rr,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
    vals = integrate.quad(vertical_integrand,-R200, R200, args=params2)
    Sigmagas[i] = rho0gas*vals[0]*(1e3*3.08e18**3/1.989e33)+1e-20
    i=i+1

#Compute Qgas
dh = hzdisc/1e4
vc = vcirc(r,dh,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)/1e5
c_s = (1.6666*1.38e-16*Tmu/1.66e-24)**0.5
kap = kappa(r,dh,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc,Tmu)
GG = 6.67e-8
Qgas = kap*c_s/np.pi/GG/(Sigmagas*1.989e33/3.08e18**2)
#print Sigmagas, Qgas, kap

#rhoa = rho_scaling(0.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)/rho_scaling(0.5,0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)
#rhob = rho_scaling(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)/rho_scaling(1.0,0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)
#rhoc = rho_scaling(2.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)/rho_scaling(2.0,0,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)
#Ta = Tmu_norm(0.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)
#Tb = Tmu_norm(1.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)
#Tc = Tmu_norm(2.0,z,rho_s,R_s,Mbulge,hbulge,Mdstar,Mdgas,hrdisc,hzdisc)

plt.figure(0)
#plt.axis([1e-3,1e2,1e-2, 1.0e4])
#plt.axis([1e-3,1e2,1e-2, 1.0e4])
plt.axis([1e-3,1e2,1e-4, 1.0e2])
#plt.axis([1e-3,1e2,1e1, 1.0e10])
plt.xscale('log')
plt.yscale('log')
plt.plot(r,vc)
#plt.plot(r,Sigmagas)
#plt.plot(r,Qgas)
#plt.plot(z,c/c*np.exp(-1))
plt.show()
