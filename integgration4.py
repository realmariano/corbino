#===================================================================================
    # Mariano Real 2019
    # Interpolation and integration using onsager for data aquired on Corbino
#===================================================================================
# import necessary Python packages
from scipy import constants
from scipy import integrate
from scipy import interpolate


import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

import os
import numpy as np
import urllib.request

# print message after packages imported successfully
print("import of packages successful")
#===================================================================================
# set working directory
#os.chdir("/home/jpalomino/earth-analytics-bootcamp/")
# print the current working directory
os.getcwd()
#===================================================================================
# Check constants and create necessary combinations of them
mu_o = constants.hbar * constants.e / (2 * constants.m_e)
wc = constants.e * 1 / (0.067 * constants.m_e)
Eo = constants.hbar * constants.e * 1 / (2 * 0.067 * constants.m_e)
print("ctte boltzman = {}".format(constants.k))
print("masa electron = {}".format(constants.m_e))
print("magneton bohr, mu_o = {}".format(mu_o))
print("cyclotron frequency, wc = {}".format(wc))
print("Eo = {}".format(Eo))

# file and workind directory, load data into datos ------------------------------------------------
archivo = "datasplineETHNov2018.txt"
datos = pd.read_csv(archivo)
# pd.read_csv gives dataFrame output, so I change it to a numpy array
p3 = datos.p3.to_numpy()
p2 = datos.p2.to_numpy()
p1 = datos.p1.to_numpy()
p0 = datos.p0.to_numpy()
B = datos.B.to_numpy()


class polinomios():
    
    def __init__(self,x0,order,*coeff):
        """ input: coeff should go in the form a_n,a_n-1,....,a1,a0
        """
        self.x0 = x0
        self.order = order
        self.coeff = list(coeff)

    def __call__(self,x):
        res = 0
        for index, coeff in enumerate(self.coeff[::-1]):
            res += coeff * (x-self.x0)**index
        return res

class dFermi():

    def __init__(self,field,temperature,alpha):
        self.field = field
        self.temperature = temperature
        self.alpha = alpha

    def __call__(self,y):
        exp1 = np.exp((y-self.field)*self.alpha/self.temperature) 
        exp2 = (np.exp((y-self.field)*self.alpha/self.temperature)+1)**2
        return exp1/exp2


# initialize arrays
xeval = []
yeval = []
idx = 0,
Binic = B[0]
# getting length of list 
length = len(B) 
temp = 0.0000001
beta = 1/(constants.k * temp)
mGaAs = 0.067
ctte_alpha = constants.hbar*constants.e * beta / (2*constants.e * mGaAs)
print(ctte_alpha)
coef_psi = ctte_alpha * beta / constants.h

for ii in range(length-1):   # Python indexes start at zero
    index = ii+1
    B1 = Binic
    B2 = B[index]
    a3=p3[index]    # because pandas import are
    a2=p2[index]
    a1=p1[index]
    a0=p0[index]
    Bx = np.linspace(B1,B2,100)

    for field in Bx:
        p = polinomios(B1,3,a3,a2,a1,a0)
        expFermi = dFermi(field,temp,ctte_alpha)
        def integrand(x):
            return expFermi(x)*p(x)*np.heaviside(x,B1)*np.heaviside(-x+B2,-B2)
        integrated = integrate.quad(integrand, B1, B2)[0] #coef_psi* integrate.quad(integrand, B1, B2)[0]
        yeval.append(integrated)
        xeval.append(field)
    
    yx = []
    Bx = []
    Binic = B2


plt.figure()
plt.plot(xeval,yeval, 'r-')
plt.legend('splines')
plt.title('Cubic-spline matlab')
plt.show()
