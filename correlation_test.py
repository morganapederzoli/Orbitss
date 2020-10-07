import numpy as np
from numpy import genfromtxt, mean
from astropy.stats import sigma_clip
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#funzioni da definire

def dist( _R, _z, _phi, _vR, _vz, _vphi, _i):
    w = np.zeros(len(_R))
    for j in range(len(_R)):
        w[j]=(sqrt((_R[j]-_R[_i])**2+(_z[j]-_z[_i])**2+(_phi[j]-_phi[_i])**2+(_vR[j]-_vR[_i])**2+(_vz[j]-_vz[_i])**2+(_vphi[j]-_vphi[_i])**2))
        #w[j]=(sqrt((_R[j]-_R[_i])**2+(_z[j]-_z[_i])**2+(_vR[j]-_vR[_i])**2+(_vz[j]-_vz[_i])**2+(_vphi[j]-_vphi[_i])**2))
    return w

#presa dati

data = genfromtxt('ComplexPlummer_full_0_E_-0.050000_Lz_0.001000.txt', delimiter='	', skip_header=1, usecols=(1,2,3,4,5,6))
R, z, phi, vR, vz, vphi = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

#superficie di sezione

f, (a,b)=subplots(2,1)

x = R[(fabs(z)<1e-4)&(vz>0)]
y = vR[(fabs(z)<1e-4)&(vz>0)]
print x.size

a.plot(x,y,'s')

step=400 #coincide con i punti considerati
count=np.zeros(step)

Max = 0
for i in range(step):
    w = dist(R,z,phi,vR,vz,vphi,i)
    if Max<= np.amax(w):
        Max=np.amax(w)
    print i
print Max, 'max'
    
hyper = np.linspace(0,int(Max)+1,step)

for i in range (x.size):

    w = dist(R,z,phi,vR,vz,vphi,i)
    hyperCopy = np.copy(hyper)

    for j in range (x.size):
        epsilon = np.amax(hyperCopy)
        for k in range(x.size):
            if w[k] < epsilon:
                count[j] += (1.0/float(len(x)**2))
        hyperCopy=np.delete(hyperCopy,-1)

rev_count=count[::-1]


b.set_yscale('log')
b.set_xscale('log')
b.plot(hyper, rev_count, 's')


plt.show()
plt.close()


#spazio dellle fasi intero

step=100 #coincide con i punti considerati
count=np.zeros(step)

Max = 0
for i in range(step):
    w = dist(R,z,phi,vR,vz,vphi,i)
    if Max<= np.amax(w):
        Max=np.amax(w)
    print i
print Max, 'max'
    
hyper = np.linspace(0,int(Max)+1,step)

for i in range(step):

    w = dist(R,z,phi,vR,vz,vphi,i)
    hyperCopy = np.copy(hyper)

    for j in range(step):

        epsilon = np.amax(hyperCopy)

        for k in range(step):

            if w[k] < epsilon:

                count[j] +=(1.0/float(step**2)) 
            
        hyperCopy=np.delete(hyperCopy,-1)

rev_count=count[::-1]

plt.yscale('log')
plt.xscale('log')

plt.plot(hyper, rev_count, 's')


plt.show()
