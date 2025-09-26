# Aharmomic potential

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

m=940
hbar=197.3
e=3.795
a=1e-2
b=10
n=1000
d=(b-a)/(n-1)

x=np.linspace(a,b,n)
c=2*m*d**2/hbar**2

#Potential
def V(x,b):
    k=100
    return 0.5*k*x**2+b*x**3/3


#Function to define matrix
def H_m(b):
    Vs=[]
    for i in range(len(x)):
        Vs.insert(i,V(x[i],b))
    H=np.empty((n,n))
    for i in range(0,n):
        H[i][i]=2+c*Vs[i]
        if i>0:
            H[i][i-1]=-1
        if i<(n-1):
            H[i][i+1]=-1
    return H

for i in [0,10,30]:
    H=H_m(i)
    Eg_val,Eg_vec=np.linalg.eigh(H)
    Energy=Eg_val[0]*hbar**2/2/m/d**2
    print('b=',i,'Energy eigen value=',Energy,'MeV')
    plt.plot(x,Eg_vec.T[0],label=(str(round(Energy,4))+"MeV"))
plt.legend()
plt.xlabel("x")
plt.grid()
plt.title("Aharmonic Potential")
plt.show()
