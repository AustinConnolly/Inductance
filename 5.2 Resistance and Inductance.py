# 5.2 Determining R and inductance

# emf vs freq

## Importing Libraries -------------------------------------------------------------------------------------------------------------------------------
import numpy as np
from scipy.optimize import curve_fit # uses L-M Algorithm
import matplotlib.pyplot as plt
from matplotlib import rc
## ---------------------------------------------------------------------------------------------------------------------------------------------------

## Reading file --------------------------------------------------------------------------------------------------------------------------------------
x = 'Induction Data.txt'
f = open(x, 'r') # opens file

n = (len(f.readlines()) - 1) # counts the number of lines in the file
f.close

f = open(x, 'r') # read and ignore header
header = f.readline()

Freqdata, Vlardata, Vaxdata = np.zeros(n), np.zeros(n), np.zeros(n) 
 
i = 0

for line in f:
    line = line.strip()
    columns = line.split()
    Freqdata[i] = float(columns[0])          # Distance from coil Data
    Vlardata[i] = 0.5*float(columns[1])        # Voltage across resistor data
    Vaxdata[i] = 0.5* (float(columns[2])*10**(-3))        # Voltage across source data
    i += 1
    
## -------------------------------------------------------------------------------------------------

## Best fit ------------------------------------------------------------------------------------------------------------
## Linearised Function -------------------------------------------------------------------------------------------------
# Data Analysis:
N = len(Freqdata)
sum_xy,sum_x,sum_y,sum_x2,sum_d2=0.0,0.0,0.0,0.0,0.0 # initialising variables

# Summation using for loops:

sum_xy = np.sum((Freqdata**2)*(Vlardata**2))
sum_x = np.sum(Freqdata**2)
sum_y = np.sum(Vlardata**2)
sum_x2 = np.sum(Freqdata**4)

# Gradient and Uncertainties:
m_num=(N*sum_xy-sum_x*sum_y) 
m_denom=(N*sum_x2-sum_x*sum_x)
m=m_num / m_denom

c_num=(sum_x2*sum_y-sum_xy*sum_x)
c_denom=(N*sum_x2-sum_x**2)
c=c_num/c_denom

for i in range(N):
        d=(Vlardata[i])**2-(m*(Freqdata[i])**2+c) # summation of d's
        sum_d2 += d**2
        
n=N/(N-2.0)

u_m = np.sqrt((sum_d2)/(N*sum_x2-sum_x**2)*n)
u_c = np.sqrt((sum_d2*sum_x2)/(N*(N*sum_x2-sum_x**2))*n)

# Prints Results:

print('m =',m,'+/-',u_m) # prints the gradient with it's uncertainty
print('c =',c,'+/-',u_c) # prints the c with it's uncertainty

print('L =',np.sqrt(m/((((np.sqrt(2))*570*10**(-3))**2)*4*(np.pi**2))))
print('R =',np.sqrt(c/((np.sqrt(2))*570*10**(-3))**2))
# Best fit
x=(Freqdata)**2
y=m*x+c

## Plotting ------------------------------------------------------------------------------------------------------------
plt.plot(Freqdata**2,Vlardata**2,'.')
plt.title('')
plt.plot(x,y,'b',label='Best fit')
#xmin,xmax,ymin,ymax =  plt.axis ([-10, 5, 0, 180])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Voltage (V)')
plt.tick_params(direction = 'in', top = 'on', right = 'on')
#plt.legend(fancybox=True,edgecolor='k',loc='upper right',framealpha=1)