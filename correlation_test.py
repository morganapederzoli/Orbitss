import numpy as np
from numpy import genfromtxt
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from fancyplot import *

#presa dati

data = genfromtxt('ComplexPlummer0_E_-0.050000_Lz_0.001000.txt', delimiter='	', skip_header=1, usecols=(1,2,3,4,5,6))
R, z, phi, vR, vz, vphi = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

n_data=R.size

#funzioni

def dist2( _R,_vR, _i,_z=np.zeros(n_data), _phi=np.zeros(n_data),  _vz=np.zeros(n_data), _vphi=np.zeros(n_data)): 
    w=(sqrt((_R-_R[_i])**2+(_z-_z[_i])**2+(_phi-_phi[_i])**2+(_vR-_vR[_i])**2+(_vz-_vz[_i])**2+(_vphi-_vphi[_i])**2))
    return w

def get_param(_x,_y):
    N=len(_x)
    #M=number of parameters
    sigma=np.zeros(N)

    #for i in range(N):
    #sigma[i]=1 	in generale (y_i-y(x_i))**/(N-M)

    sigma=np.ones(N)

    S=sum(1/sigma**2)
    S_x=sum(_x/sigma**2)
    S_y=sum(_y/sigma**2)
    S_xx=sum(_x**2/sigma**2)
    S_xy=sum(_x*_y/sigma**2)
    delta=S*S_xx - S_x**2

    q=(S_xx*S_y - S_x*S_xy)/delta
    m=(S*S_xy - S_x*S_y)/delta

    print("coeff ang :", m)
    print("normalization :", q)
    return np.array([m,q])

def min_chi_square(_x, _y):
    L=len(_x)
    l=10
    step = L - l
    param=np.zeros((step,2))
    chi_square=np.zeros(step)
    for i in range(step):
        X=_x[i:l+i]
        Y=_y[i:l+i]
        param[i]=get_param(X,Y)
        chi_square[i]=sum(Y-param[i][0]*X -param[i][1])**2 #doesn't actually work
    plt.figure()
    plt.plot(chi_square,'+')
    plt.xlabel('starting point')
    plt.ylabel('chi square')
    plt.show()
    plt.close()
    return param

#matrice distanze

distances=np.zeros((n_data,n_data)) 

for i in range(n_data):
    distances[i]= dist2(R,vR,i)

#istogramma occorrenze

distances = distances.reshape(n_data*n_data)
grid = np.logspace(-3, np.log10(np.amax(distances)), 2000)
results = np.histogram(distances, bins=grid)

#istogramma cumulativo

rr = (results[1][1:] + results[1][:-1])/2.
cumulative = np.zeros(len(rr))
cumulative[0] = results[0][0]

for i in range(1,len(rr),1):
    cumulative[i] = cumulative[i-1]+results[0][i]

#fitt

param=min_chi_square(log10(rr[100:1800]),log10(cumulative[100:1800])) 

#plottingggg

fancylayout()

f, (a,b)=subplots(2,1)

fnt = {'fontsize':16}

plt.rc('text', usetex=True)


a.set_xlabel('R', **fnt)
a.set_ylabel('vR', **fnt)
a.plot(R,vR,'s', label='Surface of Section')
a.legend(loc='upper right')

b.set_xscale('log')
b.set_yscale('log')
b.set_xlabel('r', **fnt)
b.set_ylabel('C(r)', **fnt)

b.plot(rr,cumulative,'s', label='Correlation Integral')
b.plot(rr, (10**param[500][1])*(rr**param[500][0]), label='Linear regression') #500 scelto ad occhio perchè il chi ancora non funziona
b.legend(loc='upper left')

plt.show()
