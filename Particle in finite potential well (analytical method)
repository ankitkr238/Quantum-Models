#Use the analytical solution for a particle in finite potential well

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import fsolve

l=float(input("Enter width of potential well(in Angstrom):"))
l=l*(10**-10)
V=float(input("Enter depth potential well(in eV):"))
V=V*sc.eV
m=9.1*(10**-31)
gamma=(m*V*l**2/2/(sc.hbar**2))**0.5

#using Boundary condition for even parity
def f1(x):
    '''returns eta for given xi'''
    return x*np.tan(x) 

#using Boundary condition for odd parity
def f2(x):
    '''returna eta for given xi'''
    return -x/np.tan(x)

#using Relation betwxin eta, xi and gamma
def f3(x):
    '''returns eta for given gamma and xi'''
    return np.sqrt((-(x**2)+gamma**2))


def y(x):
    return x*np.tan(x)-np.sqrt(gamma**2-x**2)

def z(x):
    return x/np.tan(x)+np.sqrt(gamma**2-x**2)

#Calculating energy eigen values 
def E(x):
    '''returns energy values for given k2'''
    return sc.hbar**2/2/m*(x**2)/sc.e

xi=[]
eta=[]

for i in range(1,4):
    p=round(fsolve(y,i)[0],4)
    q=round(fsolve(z,i)[0],4)
    if p==None and q==None:
        continue
    elif p not in xi:
        xi.append(p)
        eta.append(f1(p))
    elif q not in xi:
        xi.append(q)
        eta.append(f2(q))

xi.sort()
eta.sort(reverse=True)

#plotting transcendental equations conatining eta , xi and gamma

xxx=np.linspace(0,2*gamma,100)
xgama=np.linspace(0,gamma,100)

plt.plot(xxx,f1(xxx),color='blue',label="even parity")
plt.plot(xxx,f2(xxx),color='red',label='odd parity')
plt.plot(xgama,f3(xgama),color='green')
plt.scatter(xi,eta)
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\eta$')
plt.ylim(0,2*gamma)
plt.xlim(0,2*gamma)
plt.legend()
plt.title("Transcendental Equations")
plt.show()


#plotting eigen functions for different energy eigen values

xx1=np.linspace(-l/2,l/2,1000)
xx2=np.linspace(-l,-l/2,1000)
xx3=np.linspace(l/2,l,1000)
print("Energy Eigen Values(eV):")
for i in range(len(xi)):
    
    k1=eta[i]*2/l
    k2=xi[i]*2/l
    E_val=E(k2)
    print(E_val)
    A=np.sqrt(k1/(1+k1*l/2.0))
    if i%2==0:
        C=A*np.cos(k2*l/2.0)*np.exp(k1*l/2.0)
        plt.plot(xx1,A*np.cos(k2*xx1))
        plt.plot(xx2,C*np.e**(k1*(xx2)))
        plt.plot(xx3,C*np.e**(-k1*(xx3)))
    else:
        C=A*np.sin(k2*l/2.0)*np.exp(k1*l/2.0)
        plt.plot(xx1,A*np.sin(k2*xx1))
        plt.plot(xx2,-C*np.e**(k1*(xx2)))
        plt.plot(xx3,C*np.e**(-k1*(xx3)))

plt.axvline(x=-l/2)
plt.axvline(x=l/2)
plt.xlabel('x')
plt.ylabel('V')
plt.title("Wave function ")
plt.grid("True")
plt.show()
