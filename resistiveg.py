""" se importan los complementos"""
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from numpy.random import rand
import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
from PySpice.Spice.Netlist import Circuit
import numpy as np 
import time
#from PySpice.Unit import *
title = np.loadtxt('map_1.txt') #se ingresa el mapa a leer
mapa = np.matrix(title)
#print(mapa)
file='datos.txt'
Rc = 1
Rdia = (2*Rc)**0.5
R=0
""" se definen los puntos de inicio y fin manualmente"""
inicio = [1,1]
fin = [3,2]
V=10     
Rcont =0        
""""este es el mayor numero de pasos que puede hacer el robot para evitar que se cicle"""
dimension = mapa.shape
RMax = dimension[0]
#print(RMax)
CMax = dimension [1]
#print(CMax)
tiempo0=time.clock()
circuito = Circuit('Resistive grid')
fd = open(file,'w')
for i in range (0,RMax):
    for j in range (0,CMax):
        if (mapa[i,j]== 1):
            if (mapa[i+1,j]==1)and(i+1<=RMax):
                Rcont=Rcont+1
                R=circuito.R(Rcont,   'n'+str(i+1)+str(j+1),     'n'+str(i+2)+str(j+1), Rc)
                valores =[circuito['R'+str(Rcont)]]
                fd = open(file,'a')
                fd.write(str(valores[0])+'\n')
                
            if (mapa[i,j+1]==1)and(j+1<=CMax):
                Rcont=Rcont+1
                circuito.R(Rcont,     'n'+str(i+1)+str(j+1),      'n'+str(i+1)+str(j+2), Rc)
                valores =[circuito['R'+str(Rcont)]]
                fd = open(file,'a')
                fd.write(str(valores[0])+'\n')
                
            if (mapa[i+1,j+1]==1)and(i+1<=RMax)and(j+1<=CMax):
                Rcont=Rcont+1
                circuito.R(Rcont,     'n'+str(i+1)+str(j+1),       'n'+str(i+2)+str(j+2), Rdia)
                valores =[circuito['R'+str(Rcont)]]
                fd = open(file,'a')
                fd.write(str(valores[0])+'\n')
                
            if(i>1):
                if(mapa[i-1,j+1]==1)and(i-1>=1)and(j+1<=CMax):
                    Rcont=Rcont+1
                    circuito.R(Rcont, 'n'+str(i+1)+str(j+1),      'n'+str(i)+str(j+2), Rdia)
                    valores =[circuito['R'+str(Rcont)]]
                    fd = open(file,'a')
                    fd.write(str(valores[0])+'\n')
                                     
circuito.V('ini', 'n'+str(inicio[0])+str(inicio[1]), 'n'+str(fin[0])+str(fin[1]), V)
circuito.V('fin','n'+str(fin[0])+str(fin[1]),circuito.gnd,0)
fd.close()
simulator = circuito.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.operating_point()
fd=open('voltajes.csv','w')
for i in range (0,RMax):
    for j in range (0,CMax):
        if (mapa[i,j]== 1):
                node = analysis['n'+str(i+1)+str(j+1)]
                #print('Nodo {}: {} V'.format(str(node), float(node)))
                fd=open('voltajes.csv','a')
                fd.write('Nodo, {}, {}'.format(str(node), float(node))+'\n')
                fd.close()

