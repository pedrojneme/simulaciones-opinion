import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math


e = 2.71828182845
reflexion_der = 0.001
reflexion_izq = 0.002
omega11 = 0.05
fi11 = 0.1

tab1 = open("tabla.txt")
tab2 = open("taotabla.txt")
tab5 = open("lambda_interaccion_tabla.txt")
grado = int(3)


pob_tabla = np.loadtxt(tab1)
tao_tabla = np.loadtxt(tab2)
lambda_interaccion_tabla=np.loadtxt(tab5)

vel = np.zeros(len(pob_tabla))
vel_tao = np.zeros(len(tao_tabla))
vel_lambda_interaccion = np.zeros(len(lambda_interaccion_tabla))
poblaciones = np.zeros([len(pob_tabla),4])
tao = np.zeros([len(tao_tabla),4])
lambda_interaccion = np.zeros([len(lambda_interaccion_tabla),len(lambda_interaccion_tabla[0])-1])



omega1 = np.zeros([4,4])
fi1 = np.zeros([4,4])
sigma = np.zeros([4,4])


for i in range(0,4):
    for j in range(0,4):
        omega1[i][j]=omega11
        fi1[i][j]=fi11
        omega1[i][0]=0
        fi1[i][3]=0

for i in range(0,4):
    for j in range(0,4):
        sigma[i][j]=fi1[i][j]+fi1[j][i]+omega1[i][j]+omega1[j][i]



for i in range(0,len(pob_tabla)):
    vel[i] = math.log(pob_tabla[i][0])
    for j in range(1,len(pob_tabla[i])):
        poblaciones[i][j-1] = math.log(pob_tabla[i][j])


#son para la propaganda,
for i in range(0,len(tao_tabla)):
    vel_tao[i] = math.log(tao_tabla[i][0])

    for j in range(1,5):
        tao[i][j-1] = math.log(tao_tabla[i][j])


#son para los fi y omega.
for i in range(0,len(lambda_interaccion_tabla)):
    vel_lambda_interaccion[i] = math.log(lambda_interaccion_tabla[i][0])

    for j in range(1,len(lambda_interaccion_tabla[i])):
        lambda_interaccion[i][j-1] = math.log(lambda_interaccion_tabla[i][j])



resultados_tao = np.polyfit(vel_tao,tao,grado)

resultados_lambda_interaccion1 = np.polyfit(vel_lambda_interaccion,lambda_interaccion,grado)

resultados_lambda_interaccion = np.zeros([4,4,grado+1])




# les doy el formato que deceo:
k=0
for i in range(0,4):
    for j in range(0,4):
        for h in range(0,grado+1):
            resultados_lambda_interaccion[i][j][h] = resultados_lambda_interaccion1[h][k]
        k += 1;



def func_tao(vel,i):
    return pow(e,resultados_tao[0][i]*pow(math.log(vel),3))*pow(e,resultados_tao[1][i]*pow(math.log(vel),2))*pow(e,resultados_tao[2][i]*pow(math.log(vel),1))*pow(e,resultados_tao[3][i])


def left(vel,i):
    return (func_tao(vel,i)/(func_lambda(vel,i)+func_tao(vel,i)))*reflexion_izq

def right(vel,i):
    return (func_tao(vel,i)/(func_lambda(vel,i)+func_tao(vel,i)))*reflexion_der


def func_lambda_int(vel,i,j):
    return pow(e,resultados_lambda_interaccion[i][j][0]*pow(math.log(vel),3))*pow(e,resultados_lambda_interaccion[i][j][1]*pow(math.log(vel),2))*pow(e,resultados_lambda_interaccion[i][j][2]*pow(math.log(vel),1))*pow(e,resultados_lambda_interaccion[i][j][3])


def omega(vel,i,j):
  return 1/(func_lambda_int(vel,i,j)+func_tao(vel,j))*(omega1[i][j]/sigma[i][j])*(1-pow(e,-sigma[i][j]*func_lambda_int(vel,i,j)))

def fi(vel,i,j):
  return 1/(func_lambda_int(vel,i,j)+func_tao(vel,j))*(fi1[i][j]/sigma[i][j])*(1-pow(e,-sigma[i][j]*func_lambda_int(vel,i,j)))


def L(vel,h,P0,P1,P2,P3):
  return omega(vel,0,h)*P0 +omega(vel,1,h)*P1+omega(vel,2,h)*P2 +omega(vel,3,h)*P3

def R(vel,h,P0,P1,P2,P3):
  return fi(vel,0,h)*P0 +fi(vel,1,h)*P1 +fi(vel,2,h)*P2 +fi(vel,3,h)*P3

def P11(vel,P0,P1,P2,P3):
    return (left(vel,1)+L(vel,1,P0,P1,P2,P3))*(left(vel,2)+L(vel,2,P0,P1,P2,P3))*(left(vel,3)+L(vel,3,P0,P1,P2,P3))

def P22(vel,P0,P1,P2,P3):
    return  (right(vel,0)+R(vel,0,P0,P1,P2,P3))*(left(vel,2)+L(vel,2,P0,P1,P2,P3))*(left(vel,3)+L(vel,3,P0,P1,P2,P3))

def P33(vel,P0,P1,P2,P3):
    return (right(vel,1)+R(vel,1,P0,P1,P2,P3))*(right(vel,0)+R(vel,0,P0,P1,P2,P3))*(left(vel,3)+L(vel,3,P0,P1,P2,P3))

def P44(vel,P0,P1,P2,P3):
    return (right(vel,1)+R(vel,1,P0,P1,P2,P3))*(right(vel,2)+R(vel,2,P0,P1,P2,P3))*(right(vel,0)+R(vel,0,P0,P1,P2,P3))

def epsilon(vel,P0,P1,P2,P3):
    return P11(vel,P0,P1,P2,P3)+P22(vel,P0,P1,P2,P3)+P33(vel,P0,P1,P2,P3)+P44(vel,P0,P1,P2,P3)



y = np.arange(0.0005,1,0.01)


for i in y:
