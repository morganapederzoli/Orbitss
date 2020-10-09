import numpy as np
from numpy import genfromtxt, mean
from astropy.stats import sigma_clip
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from scipy.stats import linregress

#funzioni da definire

def dist2( _R,_vR, _i,_z=np.zeros(370), _phi=np.zeros(370),  _vz=np.zeros(370), _vphi=np.zeros(370)): 
    w=(sqrt((_R-_R[_i])**2+(_z-_z[_i])**2+(_phi-_phi[_i])**2+(_vR-_vR[_i])**2+(_vz-_vz[_i])**2+(_vphi-_vphi[_i])**2))
    return w 


#presa dati

data = genfromtxt('ComplexPlummer_full_0_E_-0.050000_Lz_0.001000.txt', delimiter='	', skip_header=1, usecols=(1,2,3,4,5,6))
R, z, phi, vR, vz, vphi = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

#superficie di sezione

f, (a,b)=subplots(2,1)

x = R[(fabs(z)<1e-4)&(vz>0)]
y = vR[(fabs(z)<1e-4)&(vz>0)]

a.plot(x,y,'s', label='Surface of Section')
a.legend(loc='upper right')

distances=np.zeros((x.size,x.size))  

distances=[dist2(x,y,i) for i in range(x.size)] 
 
hyper = np.linspace(0,int(np.amax(distances))+1,x.size) 
  
count=np.zeros(x.size) 
  
for j in range (x.size): 
    for i in range (x.size): 
         for k in range (x.size): 
             if distances[i][k] < hyper[j]: 
                 count[j] += (1.0/float(len(x)**2)) 
  
b.set_yscale('log') 
b.set_xscale('log') 
b.plot(hyper, count, 's', label='C(r) in SoS') 
b.plot(hyper,hyper) 
b.legend(loc='upper left') 
  
 
plt.show() 
plt.close()   
