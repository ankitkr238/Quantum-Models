# Coulomb Potential

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

m=0.511e6
hbar=1973
e=3.795
a=1e-2
b=10
n=1000
d=(b-a)/(n-1)

x=np.linspace(a,b,n)
c=2*m*d**2/hbar**2

#Potential
def V(x):
    return -e**2/x

Vs=[]
for i in range(len(x)):
    Vs.insert(i,V(x[i]))

H=np.empty((n,n))
for i in range(0,n):
    H[i][i]=2+c*Vs[i]
    if i>0:
        H[i][i-1]=-1
    if i<(n-1):
        H[i][i+1]=-1

Eg_val,Eg_vec=np.linalg.eigh(H)

print("Energy Eigen Values:")
for i in range(3):
    Energy=Eg_val[i]*hbar**2/2/m/d**2
    print(round(Energy,3),'eV')
    plt.plot(x,Eg_vec.T[i],label=(str(round(Energy,4))+"eV"))
plt.legend()
plt.xlabel("x")
plt.grid()
plt.title("Coulomb Potential")
plt.show()
