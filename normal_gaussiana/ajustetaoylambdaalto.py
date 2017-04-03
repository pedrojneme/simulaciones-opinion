

import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math


tab = open("taoylambda.txt","r")
tabla = np.loadtxt(tab)

vel = np.zeros(len(tabla))
tao = np.zeros(len(tabla))
landa = np.zeros(len(tabla))
e = 2.71828182845

for i in range(0,len(tabla)):
    vel[i] = tabla[i][0]
    tao[i] = tabla[i][1]
    landa[i] = tabla[i][2]

ajustex=list()
ajustey=list()

for i in range(0,len(tabla)):
    ajustex.append(math.log(vel[i]))
    ajustey.append(math.log(landa[i]))

p,q,b,n,j,k = np.polyfit(ajustex, ajustey, 5)



ajustex2=list()
ajustey2=list()
for i in range(0,len(tabla)):
    ajustex2.append(math.log(vel[i]))
    ajustey2.append(math.log(tao[i]))

a,h,c,m,l,f = np.polyfit(ajustex2, ajustey2, 5)



print "resultados cuadratico"
print "a=",a
print "h=",h
print "c=",c
print "m=",m
print "l=",l
print "f=",f
print  " "
print "p=",p
print "q=",q
print "b=",b
print "n=",n
print "j=",j
print "k=",k

def ajustelanda(vel):
    #return pow(vel,p)*pow(e,b)
    return pow(e,p*pow(math.log(vel),5))*pow(e,q*pow(math.log(vel),4))*pow(e,b*pow(math.log(vel),3))*pow(e,n*pow(math.log(vel),2))*pow(e,j*pow(math.log(vel),1))*pow(e,k)
def ajustetao(vel):
    #return pow(vel,a)*pow(e,c)
    return pow(e,a*pow(math.log(vel),5))*pow(e,h*pow(math.log(vel),4))*pow(e,c*pow(math.log(vel),3))*pow(e,m*pow(math.log(vel),2))*pow(e,l*pow(math.log(vel),1))*pow(e,f)


y = np.arange(0.0004,10,0.001)

plt.plot(vel,landa,'^',label="mean collision time")
plt.plot(y,[ajustelanda(i) for i in y],label='ajuste cuadratico')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

plt.plot(vel,tao,'s',label="mean free time" )
plt.plot(y,[ajustetao(i) for i in y],label="ajuste cuadratico")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
