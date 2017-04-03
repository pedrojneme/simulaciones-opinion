import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math


tab1 = open("tabla.txt")
pob_tabla = np.loadtxt(tab1)


plt.plot([pob_tabla[i][0] for i in range(0,len(pob_tabla))],[pob_tabla[i][1] for i in range(0,len(pob_tabla))],'.',label='PA+')
plt.plot([pob_tabla[i][0] for i in range(0,len(pob_tabla))],[pob_tabla[i][2] for i in range(0,len(pob_tabla))],'s',label='PA-',alpha=0.4 )
plt.plot([pob_tabla[i][0] for i in range(0,len(pob_tabla))],[pob_tabla[i][3] for i in range(0,len(pob_tabla))],'x',label='PB-') #B-
plt.plot([pob_tabla[i][0] for i in range(0,len(pob_tabla))],[pob_tabla[i][4] for i in range(0,len(pob_tabla))],'^',label='PB+',alpha=0.8) #B+

plt.legend(ncol=2)
plt.xlim(0.00005,2)
plt.ylim(-0.05,1.05)
plt.xscale('log')
plt.xlabel('v (arb. units)')
plt.ylabel('PA+  PA-   PB-   PB+')
plt.show()
