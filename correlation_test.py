import numpy as np
from numpy import genfromtxt, mean
from astropy.stats import sigma_clip
from matplotlib.pylab import *
import matplotlib.pyplot as plt

f, (a,b)=subplots(2,1)

data = genfromtxt('ComplexPlummer_full_0_E_-0.050000_Lz_0.001000.txt', delimiter='	', skip_header=1, usecols=(0,1,2,3,4,5,6))
t, R, z, phi, vR, vz, vphi = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:, 6]

x=[]
y=[]

for i in range(len(R)):
    if z[i]>-(10**(-4)):
        if z[i]<10**(-4):
            if vz[i]>0 :
                x.append(R[i])
                y.append(vR[i])

a.plot(x,y,'s')


def dist( _R, _z, _phi, _vR, _vz, _vphi, _i):
    w = []
    for j in range(len(_R)):
        w.append(sqrt((_R[j]-R[_i])**2+(_z[j]-z[_i])**2+(_phi[j]-phi[_i])**2+(_vR[j]-vR[_i])**2+(_vz[j]-vz[_i])**2+(_vphi[j]-vphi[_i])**2))
        #w.append(sqrt((_R[j]-R[_i])**2+(_z[j]-z[_i])**2+(_vR[j]-vR[_i])**2+(_vz[j]-vz[_i])**2+(_vphi[j]-vphi[_i])**2))

    w=np.sort(w)
    return w


step=100
count=np.zeros(step)

Max = 0
for i in range(step):
    w = dist(R,z,phi,vR,vz,vphi,i)
    for j in range(step):
        if Max <= w[j]:
            Max = w[j]
print Max
    
hyper = np.linspace(0,int(Max)+1,step)


for i in range(step):

    w = dist(R,z,phi,vR,vz,vphi,i)
    hyperCopy = np.copy(hyper)

    for j in range(step):

        epsilon = np.amax(hyperCopy)

        for k in range(step):

            if w[k] < epsilon:

                count[k] +=(1.0/float(step**2)) 
            
        hyperCopy=np.delete(hyperCopy,-1)
    

rev_count=count[::-1]
print rev_count



b.plot(hyper, rev_count, 's')


plt.show()