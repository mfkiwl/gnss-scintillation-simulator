# This script plots an intensity spectrum produced by ispectrum. It is not necessary
# to generate an intensity spectrum  spectrum.

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('ispectrum.dat')
mu_vec = data[:,0]
i_vec  = data[:,1]

params = np.loadtxt('ispectrum.log', skiprows=1)
#   Ustar         U1           U2         p1       p2      mu_break     mu_outer     mu_inner      S4        sigphi   sigNfac num

fig = plt.figure(1)
ax = fig.add_subplot(111)

plt.plot(mu_vec, i_vec, 'bo-', ms=4)
plt.rc("font", size=14)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\mu$')
plt.ylabel(r'I($\mu$)')
if params[3] == params[4]:
    title = 'U: %8.2f, p: %8.2f' \
        % (params[0],params[3],)
else:
    title = 'U: %8.2f, p$_1$: %6.2f, p$_2$: %6.2f, $\mu_b$: %8.2f' \
        % (params[0],params[3],params[4],params[5],)
if params[3] != params[4] and params[5] > 0:
    plt.axvline(x=params[5],c='r',ls='--')
if params[6] > 0:
    title = title + ', $\mu_o$: %8.2f' % (params[6],)
    plt.axvline(x=params[6],c='r',ls='--')
if params[7] > 0:
    title = title + ', $\mu_i$: %8.2f' % (params[7], )
    plt.axvline(x=params[7],c='r',ls='--')

#title = 'U: %8.2f, p$_1$: %6.2f, p$_2$: %6.2f, $\mu_b$: %8.2f, $\mu_o$: %8.2f, $\mu_i$: %8.2f' \
#        % (params[0],params[3],params[4],params[5],params[6],params[7],)
title = " ".join(title.split()) # remove duplicated spaces
plt.title(title)
label = 'S4: %8.3f' % (params[8],)
ax.annotate(label,xy=(0.9, 0.9), xycoords='axes fraction',horizontalalignment='right')

plt.show()
plt.savefig('ispectrum.png')
