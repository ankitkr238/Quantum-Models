# Using Finite difference method to solve for particle in a box and particle in harmonic potential

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

m=9.1*(10**-31)
a=-50e-10
b=50e-10
n=1000
d=(b-a)/(n)

x=np.linspace(a,b,n)
c=2*m*d**2/sc.hbar**2

#Potential
def V(x):
    k=100 #put k=0 for particle in a box
    return 0.5*k*x**2
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
    Energy=Eg_val[i]*sc.hbar**2/2/m/d**2/sc.e
    print(Energy,'eV')
    plt.plot(x,Eg_vec.T[i],label=(str(Energy)+"eV"))
plt.legend()
plt.xlabel("x")
plt.ylabel("Wave-function")
plt.axvline(x=a)
plt.axvline(x=b)
plt.title("Particle in a Box")
plt.show()

