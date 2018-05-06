""" se importan los complementos"""
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from numpy.random import rand
import PySpice.Logging.Logging as Logging
from pruebas import selector
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
""" pasamos el netlist del txt a una matriz"""
datos=np.loadtxt('datos.txt', dtype='str')
mnl1=np.matrix(datos)             #matriz del net list
mnl=mnl1.transpose()  # transpuesta 
jmax=mnl.shape[1]-1   # jmax = 6
imax=mnl.shape[0]-1   # imax = 3
n=0
m=0
ncont=0
nodos=[]
nin=[0 for x in range(jmax**2)] #nodos de entrada
nout=[0 for x in range(jmax**2)]#nodos de salida
#print (imax , jmax)
#print(mnl)
# rutina para contar nodos y eliminar repetidos
for i in range(1,imax):
    for j in range(0,jmax):
        if(mnl[i,j]!=('n'+str(fin[0])+str(fin[1]))):
            if(i<=1):
                if(mnl[i,j]!=ncont):
                    ncont= mnl[i,j]
                    nin[n]=mnl[i,j]
                    #print(mnl[i,j])
                    n=n+1   
            if(i==2):
                if(mnl[i,j]!=ncont):
                    ncont= mnl[i,j]
                    nout[m]=mnl[i,j]
                    #print(mnl[i,j])
                    m=m+1                    
for k in range (0,jmax-1):
    for l in range (0,jmax-1):
        if (nin[k] == nout[l]):
            nout.pop(l) 
      
nin.extend(nout)
nodos =[elemento for elemento in nin if elemento != 0]
nodos=set(nodos)
nodos =sorted(list(nodos))
print (nodos)
""" MNA matrix"""
mnasize = len(nodos)+1  
G=[0 for x in range (len(mnl1))]  
for i in range (0,len(mnl1)):
    G[i]=1/float(mnl[3,i])
fd=open('corrientes.csv','w')

for i in range (jmax+1):
    fd=open('corrientes.csv','a')
    if (mnl1[i,1]== 'n'+str(fin[0])+str(fin[1])):
        mnl1[i,1]= 'gnd'
        #current='I'+str(i+1),mnl1[i,1],mnl1[i,2],'G'+str(i+1)
        current='G'+str(i+1),mnl1[i,1],mnl1[i,2],str(G[i])
        
    if (mnl1[i,2]== 'n'+str(fin[0])+str(fin[1])):
        mnl1[i,2]= 'gnd'
        #current= 'I'+str(i+1),mnl1[i,1],'-'+mnl1[i,2],'G'+str(i+1)
        current= 'G'+str(i+1),mnl1[i,1],mnl1[i,2],str(G[i])
        
        
    else:
        #current= 'I'+str(i+1),mnl1[i,1],'-'+mnl1[i,2],'G'+str(i+1)
        current= 'G'+str(i+1),mnl1[i,1],mnl1[i,2],str(G[i])
        
    #print (current)
    fd.write(str(current)+'\n')
current = 'Vi','n'+str(inicio[0])+str(inicio[1]),'gnd',V
fd.write(str(current)+'\n') 
fd.close() 
  
corrientes=np.loadtxt('corrientes.csv', dtype='str')
corr=np.matrix(corrientes)
corry=corr.shape[0]
corrx=corr.shape[1]
selector ==0
if selector ==0:
    arr=[]
    row=[] 
    current=0
    for i in range(mnasize):
        row.append(0.0000000000) 
        arr.append(row)
    mna=np.matrix(arr)
    pos1=0
    pos2=0
    temp=0
    for i in range (jmax+2):
        if ('gnd' in corr[i,1])or('gnd' in corr[i,2]):
            if 'gnd' in corr[i,1]and 'Vi' not in corr[i,0] :
                for j in range (len(nodos)):
                    if nodos[j] in corr[i,2]:
                        pos2=j
                        #print('pos1 xxx')
# print('pos2 ',pos2,"  G= ", corr[i,0])
                        mna[pos2,pos2]=G[i]+mna[pos2,pos2]
            if 'gnd' in corr[i,2]and 'Vi' not in corr[i,0]:
                for j in range (len(nodos)):
                    if  nodos[j] in corr[i,1]:
                        pos1=j
#  print('pos1 ',pos1,"  G= ", corr[i,0])   
#   print('pos2 xxx')         
                        mna[pos1,pos1]=G[i]+mna[pos1,pos1] 
            if 'Vi' in corr[i,0]:
                for j in range (len(nodos)):
                    if  nodos[j] in corr[i,1]:
                        pos1=j
                        pos2=len(nodos)
                        mna[pos1,pos2]=1
                        mna[pos2,pos1]=1
                        #print (pos1,pos2,'Vi')
        else :
            for j in range (len(nodos)):
                if  nodos[j] in corr[i,1]:
                    pos1=j
                    #print('pos1 ',pos1,"  G= ", corr[i,0])
                if nodos[j] in corr[i,2]:
                    pos2=j
                    """
            if (i>0 and i<jmax):
                print(pos1,pos2,'A')
                print(mna[pos1,pos1])
                print(mna[pos2,pos2])
                print(mna[pos1,pos2])
                print(mna[pos2,pos1])
                """
                    #print('pos2 ',pos2,"  -G= ", corr[i,0])
            mna[pos1,pos1]=G[i]+mna[pos1,pos1]
            mna[pos2,pos2]=G[i]+mna[pos2,pos2]
            mna[pos1,pos2]=-G[i]+mna[pos1,pos2]
            mna[pos2,pos1]=-G[i]+mna[pos2,pos1] 
            """
            if (i>0 and i<jmax):
                print(pos1,pos2,'B')
                print(mna[pos1,pos1])
                print(mna[pos2,pos2])
                print(mna[pos1,pos2])
                print(mna[pos2,pos1])
                """
       
    fd=open('mna.csv','w')
    fd=open('mna.csv','a')
    for w in range (mnasize-1):
        fd.write(str(nodos[w])+',')
    fd.write('Iv1 ,'+'\n')
    for i in range (mnasize):
        for j in range (mnasize):
            fd.write(str(mna[i,j])+',')
        fd.write('\n')
    fd.close()
    
    print(' MATRIZ MNA NUMERICA') 
    print (mna)
    