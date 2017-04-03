

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

p,q,b = np.polyfit(ajustex, ajustey, 2)
plineal,blineal = np.polyfit(ajustex, ajustey, 1)


ajustex2=list()
ajustey2=list()
for i in range(0,len(tabla)):
    ajustex2.append(math.log(vel[i]))
    ajustey2.append(math.log(tao[i]))

a,h,c = np.polyfit(ajustex2, ajustey2, 2)
alineal,clineal = np.polyfit(ajustex2, ajustey2, 1)


print "resultados cuadratico"
print "a=",a
print "h=", h
print "c=",c
print "p=",p
print "q=", q
print "b=", b
print " "
print " "
print "resultados lineal:"
print "a=",alineal
print "c=",clineal
print "p=",plineal
print "b=", blineal


def ajustelanda(vel):
    #return pow(vel,p)*pow(e,b)
    return pow(e,p*pow(math.log(vel),2))*pow(e,q*math.log(vel))*pow(e,b)
def ajustetao(vel):
    #return pow(vel,a)*pow(e,c)
    return pow(e,a*pow(math.log(vel),2))*pow(e,h*math.log(vel))*pow(e,c)
def ajustelandalineal(vel):
    #return pow(vel,p)*pow(e,b)
    return pow(e,plineal*pow(math.log(vel),1))*pow(e,blineal)
def ajustetaolineal(vel):
    #return pow(vel,a)*pow(e,c)
    return pow(e,alineal*pow(math.log(vel),1))*pow(e,clineal)



y = np.arange(0.0001,10,0.001)

plt.plot(vel,landa,'^',label="mean collision time")
plt.plot(y,[ajustelandalineal(i) for i in y],label='ajuste lineal')
plt.plot(y,[ajustelanda(i) for i in y],label='ajuste cuadratico')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

plt.plot(vel,tao,'s',label="mean free time" )
plt.plot(y,[ajustetaolineal(i) for i in y],label="ajuste lineal")
plt.plot(y,[ajustetao(i) for i in y],label="ajuste cuadratico")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
