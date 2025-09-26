#Visualize the spherical harmonics by plotting the probability density for various values of the quantum numbers (l, m)

import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt


#Function to evaluate Angular part of Wave Function

def Y(l,m):
    
    return (-1)**m*np.sqrt(((2*l+1)*sp.factorial(l-abs(m)))/(
        4*np.pi*sp.factorial(l+abs(m))))*sp.lpmv(m,l,np.cos(theta))*np.e**(1j*m*phi)
        #used sp.lpmv to solve Assosciated legendre function

#Taking Input for Azimuthal Quantum no.

l=int(input('enter l:'))

#Meshgrid of all theta and phi values

theta=np.linspace(0,np.pi,100)
phi=np.linspace(0,2*np.pi,100)
theta,phi=np.meshgrid(theta,phi)

#Visualisation

m=-l
fig = plt.figure()

#Loop to create subplots of all Magnetic Quantum no. for a given Azimuthal Quantum no.
for i in range(2*l+1):
    Y_l_m=Y(l,m)

    #Probability Density
    P=np.abs(Y_l_m)**2

    #Creating Cartesian Coordinates
    x=P * np.sin(theta) * np.cos(phi)
    y=P * np.sin(theta) * np.sin(phi)
    z=P * np.cos(theta)

    #Plotting the points
    ax = fig.add_subplot(int(np.ceil((2*l+1)**0.5)),int(np.ceil((2*l+1)**0.5)),i+1, projection='3d')
    ax.plot_surface(x,y,z,linewidth=0,alpha=0.8,facecolors=plt.cm.viridis(P/np.max(P)), antialiased=False)
    ax.set_title(f"$|Y_{l}^{{{m}}}|^2$")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    m+=1

plt.tight_layout()
plt.show()
