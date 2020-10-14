import numpy as np
from numpy import genfromtxt, mean
from astropy.stats import sigma_clip
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from scipy.stats import linregress

#funzioni da definire

def dist2( _R,_vR, _i,_z=np.zeros(17102), _phi=np.zeros(17102),  _vz=np.zeros(17102), _vphi=np.zeros(17102)): 
    w=(sqrt((_R-_R[_i])**2+(_z-_z[_i])**2+(_phi-_phi[_i])**2+(_vR-_vR[_i])**2+(_vz-_vz[_i])**2+(_vphi-_vphi[_i])**2))
    return w

def min_chi_square(_x,_y):
    N=len(_x)
    #M=number of parameters
    sigma=np.zeros(N)

    for i in range(N):
        sigma[i]=1 # in generale (y_i-y(x_i))**/(N-M)

    S=sum(1/sigma**2)
    S_x=sum(_x/sigma**2)
    S_y=sum(_y/sigma**2)
    S_xx=sum(_x**2/sigma**2)
    S_xy=sum(_x*_y/sigma**2)
    delta=S*S_xx - S_x**2

    q=(S_xx*S_y - S_x*S_xy)/delta
    m=(S*S_xy - S_x*S_y)/delta

    print("coeff ang :", m)
    return np.array([m,q])


#presa dati

data = genfromtxt('ComplexPlummer0_E_-0.050000_Lz_0.001000.txt', delimiter='	', skip_header=1, usecols=(1,2,3,4,5,6))
R, z, phi, vR, vz, vphi = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

#superficie di sezione


distances=np.zeros((R.size,R.size)) 

for i in range(R.size):
    distances[i]= dist2(R,vR,i)

distances = distances.reshape(R.size*R.size)
grid = np.logspace(-3, np.log(np.amax(distances)), 200)
results = np.histogram(distances, bins=grid)

rr = (results[1][1:] + results[1][:-1])/2.
cumulative = np.zeros(len(rr))
cumulative[0] = results[0][0]

for i in range(1,len(rr),1):
    cumulative[i] = cumulative[i-1]+results[0][i]

param=min_chi_square(log(rr[:130]),log(cumulative[:130])) #point 130 should exclude final plateu



f, (a,b)=subplots(2,1)

#x = R[(fabs(z)<1e-4)&(vz>0)]
#y = vR[(fabs(z)<1e-4)&(vz>0)]


fnt = {'fontsize':16}

plt.rc('text', usetex=True)

a.plot(R,vR,'s', label='Surface of Section')
a.legend(loc='upper right')

b.set_xscale('log')
b.set_yscale('log')
b.plot(rr,cumulative,'s')
b.plot(rr, (10**param[1])*(rr**param[0]))

b.set_xlabel('R $\\rho$', **fnt)

plt.show() #pirla se non lo metti non mostra svegliati
plt.close()
