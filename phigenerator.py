import matplotlib.pyplot as plt
import scipy.integrate as scint
import scipy.interpolate as scintpol
import scipy.special as scsp
import numpy as np

alpha, Te = np.loadtxt('inputfile.txt')
print(alpha, Te)
Ts = 1.0 + Te

def xifp(x, phi):
	y = 1.0/np.sqrt( -3.0 -2.0*phi + 4.0*np.exp(phi) - np.exp(2.0*phi) )
	return y	

xarr2 = np.arange(1.0, 10.0, 0.05)
garr1 = np.arange(0.0, 1.0, 0.05)
xarr1 = [i**2 for i in garr1]
xarr = np.concatenate((xarr1, xarr2))
phiarr = np.arange(np.log(alpha), -0.0001, 0.0001)
phiarrf = [p + 0.5*alpha**2*(np.exp(-2.0*p)-1.0) for p in phiarr]
xif = scint.odeint(xifp, 0.0, phiarrf)
xif  = np.resize(xif, len(xif))
f1 = scintpol.interp1d(xif, phiarr)
phi1 = f1(xarr)

fp = open('phidata.txt', 'w')

xtop = 15.0*(np.sqrt(Ts))
xcr=1.0*np.sqrt(Ts)
deltax = 0.08*np.sqrt(Ts)
nn1 = (int) ( 2.0*np.sqrt(xcr)/deltax )
nn2 = (int) ( (xtop - xcr)/deltax )
for i in range(nn1):
	g = i*np.sqrt(xcr)/nn1 
	fp.write(str(g)+' '+str(Te*f1(g**2/np.sqrt(Te/2)))+'\n')
	#fp.write(str(g)+' '+str(0.0)+'\n')

for i in range(nn2):
	g = np.sqrt(xcr + i*deltax)
	fp.write(str(g)+' '+str(Te*f1(g**2/np.sqrt(Te/2)))+'\n')
	#fp.write(str(g)+' '+str(0.0)+'\n')

fp.close()
