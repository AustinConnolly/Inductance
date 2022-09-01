# Plotting amplitude of emf vs distance from the centre of the coil
# CNNAUA001
# 4/10/2020

## Importing Libraries -------------------------------------------------------------------------------------------------------------------------------
import numpy as np
from scipy.optimize import curve_fit # uses L-M Algorithm
import matplotlib.pyplot as plt
from matplotlib import rc
## ---------------------------------------------------------------------------------------------------------------------------------------------------

## Reading file --------------------------------------------------------------------------------------------------------------------------------------
x = 'field Axis Data.txt'
f = open(x, 'r') # opens file

n = (len(f.readlines()) - 1) # counts the number of lines in the file
f.close

f = open(x, 'r') # read and ignore header
header = f.readline()

Distdata, Vdata = np.zeros(n), np.zeros(n) 
 
i = 0

for line in f:
    line = line.strip()
    columns = line.split()
    Distdata[i] = (float(columns[0])*10**(-2))          # Distance from coil Data
    Vdata[i] = 0.5*(float(columns[1])*10**(-3))         # Voltage across resistor data
    #SorcVdata[i] = float(columns[3])        # Voltage across source data
    i += 1
    
f.close()
## ---------------------------------------------------------------------------------------------------------------------------------------------------

## Defining function ---------------------------------------------------------------------------------------------------------------------------------
def B1(V): # prop constant
    omega = 2000*np.pi
    N = 175
    r = ((1.3*10**(-2))/2)
    b = V/((N)*(np.pi*r**2)*omega)
    return b
def B2(z): # Theoretical
    mu0 = 4*np.pi*10**(-7)
    I0 = 500*10**(-3)
    N = 120
    a = (0.5*(6.8*10**(-2)))
    b1 = ((mu0*(N)*(I0))/2)
    b2 = (a**2/((a**2+z**2)**(3/2)))
    b = b1*b2
    return b
def emf(z): # Theoretical
    mu0 = 4*np.pi*10**(-7)
    I0 = 500*10**(-3)
    N1 = 175
    N2 = 120
    r = ((1.3*10**(-2))/2)
    A = np.pi*r**2
    a = (0.5*(6.8*10**(-2)))
    omega = 2000*np.pi
    b1 = (omega*mu0*N1*N2*I0*A)/2
    b2 = (a**2/((a**2+z**2)**(3/2)))
    b = b1*b2
    return b
## ---------------------------------------------------------------------------------------------------------------------------------------------------

## Plotting graphs -----------------------------------------------------------------------------------------------------------------------------------
## emf:
#plt.plot(Distdata,Vdata,'r',label = 'emf')
#plt.title('Amplitude of emf of coil vs Distance from coil')
#xmin,xmax,ymin,ymax =  plt.axis ([-0.10, 0.05, 0, 0.16])
#plt.xlabel('Distance (m)')
#plt.ylabel('Voltage (V)')
## B fields:
plt.plot(Distdata,B1(Vdata),'b',label = 'Experimental')
plt.plot(Distdata,0.9611112*B2(Distdata),'g',label = 'Theoretical')
xmin,xmax,ymin,ymax =  plt.axis ([-0.10, 0.05, 0, 0.0011])
plt.xlabel('Distance (m)')
plt.ylabel('Magnetic Field (T)')
plt.title('Experimental Magnetic Field vs Theoretical Magnetic Field')
## emf theoretical vs experimental:
#plt.plot(Distdata, 0.9611112*emf(Distdata),'b', label = 'Theoretical')
#plt.plot(Distdata,Vdata,'g',label = 'Experimental')
#xmin,xmax,ymin,ymax =  plt.axis ([-0.10, 0.05, 0, 0.16])
#plt.xlabel('Distance (m)')
#plt.ylabel('emf (V)')
#plt.title('Experimental emf vs Theoretical emf')

plt.tick_params(direction = 'in', top = 'on', right = 'on')
plt.legend(fancybox=True,edgecolor='k',loc='upper left',framealpha=1)
## ---------------------------------------------------------------------------------------------------------------------------------------------------