#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
df=np.loadtxt("./gs_check.dat")

fig=plt.figure(figsize=(8,6))
ax=fig.add_subplot(111)
ax.plot(df[:,0],df[:,1],marker='+',label='Lanczos',ls='none',c='r',mew=1,ms=10)
ax.plot(df[:,0],df[:,2],marker='o',label='Lieb-Wu',ls='none',c='b',mew=1,ms=10,mfc='none')
ax.set_xlabel(r'$U$',fontsize=24)
ax.set_xlim(-0.5,61)
ax.set_ylabel(r'$E_{gs}$',fontsize=24)

plt.legend(loc='lower right',numpoints=1,fontsize=20)
plt.show()
