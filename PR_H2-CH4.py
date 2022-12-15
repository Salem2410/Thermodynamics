import numpy as np
print("################# Number of Molecular Units in a Cubic Simulation Box using Peng-Robinson EoS for H2-CH4 Mixture #################")
#Prompt for mole fractions
T = float(input('Temperature in Kelvin : ')) 
P = float(input('Pressure in psi : '))
P = P *  6894.7572931783 #Pa
al =  float(input('Simulation box side length (nm) : '))      
V = (al*10**-9)**3 #m3 
yH2 = float(input('Mole fraction of H2 : '))
if yH2 > 1:
    quit("Mole fraction must be less than or equal 1")
yCH4 = 1 - yH2
#H2 properties
Tc_H2 = 33.19 #K 
Pc_H2 = 13.13*10**5#Pa
w_H2 = -0.216
#CH4 properties
Tc_CH4 = 190.6 #K 
Pc_CH4 = 45.99*10**5 #Pa
w_CH4 = 0.012
#Constants
R = 8.314 #J/mol/K
NA = 6.02214076e+23 #mol^-1
sigma = 1 + np.sqrt(2)
epsilon = 1 - np.sqrt(2)
Omega = 0.07780
Psi = 0.45724
#Reduced temperature
Tr_H2 = T/Tc_H2
Tr_CH4 = T/Tc_CH4
#Define alpha function
def alpha(Tr, w):
    return (1+(0.37464+1.54226*w-0.26992*w**2)*(1-np.sqrt(Tr)))**2
#a
a_H2 = Psi*alpha(Tr_H2, w_H2)*R**2*Tc_H2**2/Pc_H2
a_CH4 = Psi*alpha(Tr_CH4, w_CH4)*R**2*Tc_CH4**2/Pc_CH4
#b
b_H2 = Omega*R*Tc_H2/Pc_H2
b_CH4 = Omega*R*Tc_CH4/Pc_CH4
#Combining rule
aij = np.sqrt(a_H2*a_CH4)
#Mixing rules
a = yH2**2*a_H2+2*yH2*yCH4*aij+yCH4**2*a_CH4
b = yH2*b_H2+yCH4*b_CH4
#Numerical solution
Beta = b*P/R/T
q = a/b/R/T
Z = 1
tol = 1
while tol > 0.0001:
    t1 = 1 + Beta    
    t2 = q*Beta*(Z-Beta)/(Z+sigma*Beta)/(Z+epsilon*Beta)
    Zn = t1 - t2
    tol = abs(Zn-Z)
    Z = Zn
v = Z*R*T/P #m3/mol
n = V/v
N = n*NA
print("Number of H2 molecules:",round( N*yH2,1))
print("Number of CH4 molecules:",round(N*yCH4,1))
