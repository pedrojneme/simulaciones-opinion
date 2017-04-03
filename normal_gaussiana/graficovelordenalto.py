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

r1 = 0.001
l1 = 0.002
w1 = 0.05
f1 = 0.1
sigma = 2*w1 + 2*f1
e = 2.71828182845



a= -2.97524351468e-05
h= -6.87592732823e-05
c= -0.00388213952617
m= -0.00430405003045
l= -0.987006232059
f= 2.40796548125

p= -3.19710083334e-07
q= 0.00197536072825
b= 0.0209104237593
n= 0.0033731799371
j= -1.0049348705
k= 0.0668191592609








def tao(vel):
  return pow(e,a*pow(math.log(vel),5))*pow(e,h*pow(math.log(vel),4))*pow(e,c*pow(math.log(vel),3))*pow(e,m*pow(math.log(vel),2))*pow(e,l*pow(math.log(vel),1))*pow(e,f)

def landa(vel):
  return pow(e,p*pow(math.log(vel),5))*pow(e,q*pow(math.log(vel),4))*pow(e,b*pow(math.log(vel),3))*pow(e,n*pow(math.log(vel),2))*pow(e,j*pow(math.log(vel),1))*pow(e,k)

def omega(vel):
    return (1/(tao(vel)+landa(vel)))*((w1/sigma)*(1-pow(e,-sigma*landa(vel))))
def fi(vel):
    return (1/(tao(vel)+landa(vel)))*((f1/sigma)*(1-pow(e,-sigma*landa(vel))))
def left(vel):
    return (tao(vel)/(tao(vel)+landa(vel)))*l1
def right(vel):
    return  (tao(vel)/(tao(vel)+landa(vel)))*r1
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

for i in range(len(tabla)):
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
