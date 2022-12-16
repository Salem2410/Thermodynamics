import numpy as np
natural = tuple('123456789')
n = input("Number of components: ")
EoS = str(input("Equation of state (vdW, RK, SRK, or PR): "))
if EoS != "vdW" and EoS != "RK" and EoS != "SRK" and EoS != "PR":
    quit("Invalid entry for equation state, enter 'vdW', 'RK', 'SRK' or 'PR' instead")
if not(n in natural):
    quit("Invalid input: number of components must be natural (e.g. 1, 2, etc...)")
n = int(n)
T = float(input("Temperature (in Kelvin): "))
P = float(input("Pressure (in bar): ")) * 10**5
R = 8.314 #J/mol/K
Comp = []
Tc = []
Pc = []
w = []
y = []
a = []
amix = []
b = []
bmix = []
def ab(EoS, T, Tc, Pc, w):
    Tr = T/Tc #Reduced temperature
    R = 8.314 #J/mol/K
    if EoS == 'vdW':
        Omega = 1/8
        Psi = 27/64
        alpha = 1
    elif EoS == 'RK':
        Omega = 0.08664
        Psi = 0.42748
        alpha = np.sqrt(Tr)**-1
    elif EoS == 'SRK':
        Omega = 0.08664
        Psi = 0.42748
        alpha = (1+(0.480+1.574*w-0.176*w**2)*(1-np.sqrt(Tr)))**2
    elif EoS == 'PR':
        Omega = 0.07780
        Psi = 0.45724
        alpha = (1+(0.37464+1.54226*w-0.26992*w**2)*(1-np.sqrt(Tr)))**2
    a = Psi*alpha*R**2*Tc**2/Pc
    b = Omega*R*Tc/Pc
    return a, b
for i in range(0,n):
    comp_name = str(input("Name of component %i: " % (i+1))) 
    if i != n-1:
        yi = float(input("Mole fraction for "+comp_name+": "))
    else:
        if sum(y) >=1:
            quit("Wrong mole fractions")
        yi = 1 - sum(y)
    Tci = float(input("Critical temperature (in Kelvin) for "+comp_name+": "))
    Pci = float(input("Critical pressure (in bar) for "+comp_name+": ")) * 10**5
    wi = float(input("Acentric factor for "+comp_name+": ")) 
    ai, bi = ab(EoS, T, Tci, Pci, wi)
    Comp.append(comp_name)
    Tc.append(Tci)
    Pc.append(Pci)
    w.append(wi)
    y.append(yi)
    a.append(ai)
    b.append(bi)
for i in range(0,n):
    for j in range(0,n):
        amix_i = y[i]*y[j]*np.sqrt(a[i]*a[j])
        bmix_i = y[i]*y[j]*np.sqrt(b[i]*b[j])
        amix.append(amix_i)
        bmix.append(bmix_i)
am = sum(amix)
bm = sum(bmix)
if EoS == 'vdW':
    sigma = 0
    epsilon = 0
elif EoS == 'RK':
    sigma = 1 
    epsilon = 0
elif EoS == 'SRK':
    sigma = 1 
    epsilon = 0
elif EoS == 'PR':
    sigma = 1 + np.sqrt(2)
    epsilon = 1 - np.sqrt(2)
Beta = bm*P/R/T
q = am/bm/R/T
Z = 1
tol = 1
while tol > 0.01:
    Zn = 1 + Beta - q*Beta*(Z-Beta)/(Z+sigma*Beta)/(Z+epsilon*Beta)
    tol = abs(Zn-Z)/Zn*100
    Z = Zn
V = Z*R*T/P
print("Z =", round(Z,3))
print("V =", round(V,3), "m3/mol =", round(V*10**6,3), "cm3/mol")
