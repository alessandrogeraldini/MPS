import numpy as np
import scipy.integrate as int
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
import scipy.special as scsp

Tevect = [x*0.01 for x in range(1, 99)]
for i in range(6):
	print(i)
	Tevect.append(1 - (0.1-0.02*i)**2)

print(Tevect)
kvect = []

for i in range(len(Tevect)):
	Te = Tevect[i]
	print(Te)
	k = 0.0
	kold = 1000.0
	while abs(k - kold) > 0.000001:
		kold = k
		lhs = 4*(np.exp(k**2/2)*scsp.i0(k**2/2) -1)
		rhs = 2/Te - 2
		print(k, lhs, rhs)
		#result = lhs - rhs
		k = np.sqrt( 2.0*np.log( ( rhs/4 + 1 )/scsp.i0(kold**2/2) ) )
	
	kvect.append(k)
	#kvect.append(float(format(k, '.4f')))


print(len(Tevect), len(kvect))
print(Tevect, kvect)
Tivect = [1/x for x in Tevect]

plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Palatino']})
plt.xlabel(r'$\tau$', fontsize=26)
plt.ylabel(r'$ k \rho_{i} $', fontsize=26)
plt.xlim(0.0, 20.0)
plt.ylim(0.0, 2.0)
###plt.title(r'$T_e / T_i = 1$')
plt.xticks(range(21))
plt.tight_layout()
np.savetxt('k1gyrokinetic.txt', zip(Tevect, kvect), fmt='%.5f')
plt.plot(Tivect, kvect)

plt.savefig('fig-PTD-k1vstau.pdf', dpi =2000)
