#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import cauchy

fig, ax = plt.subplots(1,1)
mean,var,skew,kurt = cauchy.stats(moments='mvsk')

x = np.linspace(cauchy.ppf(0.01),cauchy.ppf(0.99),1000)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax.plot(x, cauchy.pdf(x), color='#2c0032', lw=3, label='Standardized')
ax.plot(x, cauchy.pdf(x,loc=0,scale=0.5), color='#8c07ee', lw=3, label='$\Delta$ = 0.5')
#ax.plot(x, cauchy.pdf(x), color='#2f0059', lw=3, label='Standardized')
ax.plot(x, cauchy.pdf(x,loc=0,scale=2.0), color='#ad1927', lw=3, label='$\Delta$ = 2.0')
ax.plot(x, cauchy.pdf(x,loc=0,scale=3.0), color='#f5ce00', lw=3, label='$\Delta$ = 3.0')

ax.legend(loc='best',frameon=False)
ax.set_xlim([-13.0,13.0])
#ax.set_title('Cauchy PDF')

plt.text(-11,0.5,'$f(z) = \\frac{\Delta}{\pi(z^2 + \Delta^2)}$', fontsize=20)
fig.savefig('cauchy.eps',format='eps',dpi=1000)
plt.show()
