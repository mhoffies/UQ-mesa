#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import cauchy

fig, ax = plt.subplots(1,1)
mean,var,skew,kurt = cauchy.stats(moments='mvsk')

x = np.linspace(cauchy.ppf(0.01),cauchy.ppf(0.99),1000)
y = np.linspace(0,40,100)
ons = np.ones(100)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax.plot(x, cauchy.pdf(x, loc=0.53432, scale=0.0382 ), color='#2c0032', lw=3, label='$\Delta=0.0382$')
ax.plot(x, cauchy.pdf(x, loc=0.55334, scale=0.0102 ), color='#8c07ee', lw=3, label='$\Delta=0.0102$')
ax.plot(x, cauchy.pdf(x, loc=0.54545, scale=0.0024 ), color='#ad1927', lw=3, label='$\Delta=0.0024$')
ax.plot(x, cauchy.pdf(x, loc=0.59613, scale=0.0988 ), color='#f5ce00', lw=3, label='$\Delta=0.0988$')
ax.plot(x, cauchy.pdf(x, loc=0.53622, scale=0.0349 ), color='#ff6600', lw=3, label='$\Delta=0.0349$')
ax.plot(x, cauchy.pdf(x, loc=0.5237, scale=0.0295 ), color='#6ca6cd', lw=3, label='$\Delta=0.0295$')
ax.plot(0.695*ons,y,'k--',lw=1)
ax.plot(0.494*ons,y,'k--',lw=1)
ax.legend(loc='best',frameon=False)
ax.set_xlim([0.35,0.85])

#ax.set_title('Cauchy PDF')

#plt.text(-2.0,0.2,'$\Delta=0.53$',fontsize=12)
#plt.text(-0.4,0.27,'$\Delta=0.55$',fontsize=12)
#plt.text(0.5,0.3,'$\Delta=0.54$',fontsize=12)
#plt.text(1.7,0.25,'$\Delta=0.59$',fontsize=12)

#plt.text(-11,0.0,'$f(z) = \\frac{\Delta}{\pi(z^2 + \Delta^2)}$', fontsize=20)

fig.savefig('manycauchy.eps',format='eps',dpi=1000)
plt.show()
