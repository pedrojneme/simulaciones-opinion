import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math

tab = open("tabla.txt","r")
tabla = np.loadtxt(tab)

x = np.zeros(len(tabla))
a1 = np.zeros(len(tabla))
b1 = np.zeros(len(tabla))
c1 = np.zeros(len(tabla))
d1 = np.zeros(len(tabla))

r = 0.001
l = 0.002
w = 0.05
f = 0.1
sigma = 2*w + 2*f
e = 2.71828182845


a= -1.00220012756
c= 2.40517223872
p= -0.926173288973
b= 0.0634758837348









def tao(vel):
    return pow(e,a*pow(math.log(vel),1))*pow(e,c)

def landa(vel):
    return pow(e,p*pow(math.log(vel),1))*pow(e,b)
def omega(vel):
    return (1/(tao(vel)+landa(vel)))*((w/sigma)*(1-pow(e,-sigma*landa(vel))))
def fi(vel):
    return (1/(tao(vel)+landa(vel)))*((f/sigma)*(1-pow(e,-sigma*landa(vel))))
def left(vel):
    return (tao(vel)/(tao(vel)+landa(vel)))*l
def right(vel):
    return  (tao(vel)/(tao(vel)+landa(vel)))*r
def epsilon(vel):
    return pow((omega(vel)+left(vel)),3)+pow((omega(vel)+left(vel)),2)*pow((fi(vel)+right(vel)),1)+pow((omega(vel)+left(vel)),1)*pow((fi(vel)+right(vel)),2)+pow((fi(vel)+right(vel)),3)


def PAA(vel):
    return (1/epsilon(vel))*pow((omega(vel)+left(vel)),3)

def PA(vel):
    return (1/epsilon(vel))*pow((omega(vel)+left(vel)),2)*pow((fi(vel)+right(vel)),1)

def PB(vel):
    return (1/epsilon(vel))*pow((omega(vel)+left(vel)),1)*pow((fi(vel)+right(vel)),2)

def PBB(vel):
    return (1/epsilon(vel))*pow((fi(vel)+right(vel)),3)

y = np.arange(0.0001,10,0.001)

for i in range(0,len(tabla)):
 x[i] = tabla[i][0]
 a1[i] = tabla[i][1]
 b1[i] = tabla[i][2]
 c1[i] = tabla[i][3]
 d1[i] = tabla[i][4]

#grafico:
plt.plot(y,[PAA(i) for i in y],label='PA+ ')
plt.plot(y,[PA(i) for i in y],label='PA- ')
plt.plot(y,[PB(i) for i in y],label='PB- ')
plt.plot(y,[PBB(i) for i in y],label='PB+ ')
plt.plot(x,a1,'.',label='PA+')
plt.plot(x,b1,'s',label='PA-',alpha=0.4 )
plt.plot(x,c1,'x',label='PB-') #B-
plt.plot(x,d1,'^',label='PB+',alpha=0.8) #B+

plt.legend(ncol=2)
plt.xlim(0.0004,10)
plt.ylim(-0.05,1.05)
plt.xscale('log')
plt.xlabel('mean vel (arb. units)')
plt.ylabel('PA+  PA-   PB-   PB+')
plt.show()
