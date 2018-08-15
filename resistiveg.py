""" se importan complementos"""
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy import  misc
import matplotlib.pyplot as plt
import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
from PySpice.Spice.Netlist import Circuit
import numpy as np 
import time
""" #se ingresa el mapa a leer========================================================="""
title = np.loadtxt('map_1.txt') 
mapa = np.matrix(title)
file='datos.txt'
Rc = 1
Rdia = (2*Rc)**0.5
Rt=1
""" se definen los puntos de inicio y fin manualmente================================="""
inicio = [1,1]
fin = [9,9]
V=10     
Rcont =0        
""""este es el mayor numero de pasos que puede hacer el robot para evitar que se cicle"""
dimension = mapa.shape
RMax = dimension[0]
CMax = dimension [1]
""" inicia el cronometro"""
tiempo0=time.clock()
""" creacion del circuito partiendo del mapa=========================================="""
circuito = Circuit('Resistive grid')
fd = open(file,'w')
for i in range (0,RMax):
    for j in range (0,CMax):
        if (mapa[i,j]== 1):
            if (mapa[i+1,j]==1)and(i+1<=RMax):
                Rcont=Rcont+1
                R=circuito.R(Rcont,   'n'+str(i+1)+'a'+str(j+1),     'n'+str(i+2)+'a'+str(j+1), Rc)
                valores =[circuito['R'+str(Rcont)]]
                fd = open(file,'a')
                fd.write(str(valores[0])+'\n')
                
            if (mapa[i,j+1]==1)and(j+1<=CMax):
                Rcont=Rcont+1
                circuito.R(Rcont,     'n'+str(i+1)+'a'+str(j+1),      'n'+str(i+1)+'a'+str(j+2), Rc)
                valores =[circuito['R'+str(Rcont)]]
                fd = open(file,'a')
                fd.write(str(valores[0])+'\n')
                
            if (mapa[i+1,j+1]==1)and(i+1<=RMax)and(j+1<=CMax):
                Rcont=Rcont+1
                circuito.R(Rcont,     'n'+str(i+1)+'a'+str(j+1),       'n'+str(i+2)+'a'+str(j+2), Rdia)
                valores =[circuito['R'+str(Rcont)]]
                fd = open(file,'a')
                fd.write(str(valores[0])+'\n')
                
            if(i>1):
                if(mapa[i-1,j+1]==1)and(i-1>=1)and(j+1<=CMax):
                    Rcont=Rcont+1
                    circuito.R(Rcont, 'n'+str(i+1)+'a'+str(j+1),      'n'+str(i)+'a'+str(j+2), Rdia)
                    valores =[circuito['R'+str(Rcont)]]
                    fd = open(file,'a')
                    fd.write(str(valores[0])+'\n')
                                     
circuito.V('ini', 'n'+str(inicio[0])+'a'+str(inicio[1]), 'n'+str(fin[0])+'a'+str(fin[1]), V)
circuito.V('fin','n'+str(fin[0])+'a'+str(fin[1]),circuito.gnd,0)
fd.close()
simulator = circuito.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.operating_point()
fd=open('voltajes.csv','w')
for i in range (0,RMax):
    for j in range (0,CMax):
        if (mapa[i,j]== 1):
                node = analysis['n'+str(i+1)+'a'+str(j+1)]
                fd=open('voltajes.csv','a')
                fd.write('Nodo, {}, {}'.format(str(node), float(node))+'\n')
                fd.close()


""" transferencia  del netlist del txt a una matriz=========================================="""
datos=np.loadtxt('datos.txt', dtype='str')
mnl1=np.matrix(datos)
mnl=mnl1.transpose() 
jmax=mnl.shape[1]-1   
imax=mnl.shape[0]-1   
n=0
m=0
ncont=0
nodos=[]
nin=[0 for x in range(jmax**2)] 
nout=[0 for x in range(jmax**2)]

"""  contar nodos y eliminar repetidos======================================================="""
for i in range(1,imax):
    for j in range(0,jmax):
        if(mnl[i,j]!=('n'+str(fin[0])+'a'+str(fin[1]))):
            if(i<=1):
                if(mnl[i,j]!=ncont):
                    ncont= mnl[i,j]
                    nin[n]=mnl[i,j]
                    
                    n=n+1   
            if(i==2):
                if(mnl[i,j]!=ncont):
                    ncont= mnl[i,j]
                    nout[m]=mnl[i,j]
                    
                    m=m+1                    
for k in range (0,jmax-1):
    for l in range (0,jmax-1):
        if (nin[k] == nout[l]):
            nout.pop(l) 
      
nin.extend(nout)
nodos =[elemento for elemento in nin if elemento != 0]
nodos=set(nodos)
nodos =sorted(list(nodos))

""" MNA matrix======================================================================="""
mnasize = len(nodos)+1  
G=[0 for x in range (len(mnl1))]  
for i in range (0,len(mnl1)):
    G[i]=1/float(mnl[3,i])
fd=open('corrientes.csv','w')

for i in range (jmax+1):
    fd=open('corrientes.csv','a')
    if (mnl1[i,1]== 'n'+str(fin[0])+'a'+str(fin[1])):
        mnl1[i,1]= 'gnd'
        current='G'+str(i+1),mnl1[i,1],mnl1[i,2],str(G[i])
        
    if (mnl1[i,2]== 'n'+str(fin[0])+'a'+str(fin[1])):
        mnl1[i,2]= 'gnd'
        current= 'G'+str(i+1),mnl1[i,1],mnl1[i,2],str(G[i])
        
        
    else:
        current= 'G'+str(i+1),mnl1[i,1],mnl1[i,2],str(G[i])

    fd.write(str(current)+'\n')
current = 'Vi','n'+str(inicio[0])+'a'+str(inicio[1]),'gnd',V
fd.write(str(current)+'\n') 
fd.close() 
  
corrientes=np.loadtxt('corrientes.csv', dtype='str')
corr=np.matrix(corrientes)
corry=corr.shape[0]
corrx=corr.shape[1]
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

                    mna[pos2,pos2]=G[i]+mna[pos2,pos2]
        if 'gnd' in corr[i,2]and 'Vi' not in corr[i,0]:
            for j in range (len(nodos)):
                if  nodos[j] in corr[i,1]:
                    pos1=j
         
                    mna[pos1,pos1]=G[i]+mna[pos1,pos1] 
        if 'Vi' in corr[i,0]:
            for j in range (len(nodos)):
                if  nodos[j] in corr[i,1]:
                    pos1=j
                    pos2=len(nodos)
                    mna[pos1,pos2]=1
                    mna[pos2,pos1]=1
    else :
        for j in range (len(nodos)):
            if  nodos[j] in corr[i,1]:
                pos1=j

            if nodos[j] in corr[i,2]:
                pos2=j

        mna[pos1,pos1]=G[i]+mna[pos1,pos1]
        mna[pos2,pos2]=G[i]+mna[pos2,pos2]
        mna[pos1,pos2]=-G[i]+mna[pos1,pos2]
        mna[pos2,pos1]=-G[i]+mna[pos2,pos1] 
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

"""seleccionador matrices normal vs dispersa==============================================="""
selectsolve = 1
if selectsolve == 0:
    MNA=np.matrix(mna)
    mov=mnasize-1
    roff=[] 
    for i in range(mnasize):
        roff.append(0) 
    b=np.matrix(roff)
    b=np.transpose(b)
    for i in range (mnasize):
        b[mov-i,0]=(mna[i,mov])*V
    print(b)
    result=(MNA.I)*b
    print( 'resultados \n',result)
else :
    mov=mnasize-1
    roff=[] 
    for i in range(mnasize):
        roff.append(0) 
    b=np.matrix(roff)
    b=np.transpose(b)
    for i in range (mnasize):
        b[mov-i,0]=(mna[i,mov])*V
    """ resolviendo con sparce"""
    A=lil_matrix(mna)
    A = A.tocsr()
    x = spsolve(A,b)
    
    fd=open('vsp.txt','w')
    fd=open('vsp.txt','a')
    print("Lista de voltajes : \n")
    for i in range (mnasize):
        if (i < mov):
            """ guardando los datos de voltajes en un archivo ===================================="""
            fd.write(str(nodos[i])+' '+str(x[i])+'\n')
            print (nodos[i], '= V',x[i])
        else:
            print ('Iv1 = A',x[i] )
    fd.close()
""" pasando los datos del archivo a una matriz para su manipulacion"""
datax=np.loadtxt('vsp.txt', dtype='str')
""" llenando la matriz Y con los voltajes segun su nodo==========================================="""
NVoltajes=x.shape
arr=[]
row=[] 
for i in range(CMax):
    row.append(0.0)
    arr.append(row)
Y=np.matrix(arr)
print ("numero de voltajes ",NVoltajes)
mvs=np.matrix(datax)             
for i in range (mvs.shape[0]):
    xaux = mvs[i,0].split('n')
    aux = str(xaux)
    b=aux.split('a')
    bob=str(b)
    wui = bob[8:9]
    yei = bob[13:14]
    j=int(wui)-1
    k=int(yei)-1
    Y[j,k]=mvs[i,1]
print("mapa de voltajes \n",Y)   

""" LCC======================================================================================"""
xaxis=[]
yaxis=[]
i=inicio[0]-1
j=inicio[1]-1
i1=inicio[0]-1
j1=inicio[1]-1
temp=[inicio[0]-1,inicio[1]-1]
val=-100
ruta=[0,0]
temp2=ruta
newnode=[0,0]
NVoltajes=x.shape
for cont in range(NVoltajes[0]):
    i=i1
    j=j1
    nods= [[i-1, j], [i-1, j+1], [i, j+1], [i+1, j+1], [i+1, j], [i+1, j-1], [i, j-1],[i-1, j-1]]
    for k in range (7):
        nod = nods[k]
        if (nod[0]>=0 or nod[1]>=0):
            if(i!=fin[0])and(j!=fin[1]):
                if not ((nod[0]==temp[0])and(nod[1]==temp[1])):
                    if(nod[0]<= RMax) and (nod[1]<= CMax):
                        if (mapa[nod[0],nod[1]]==1):
                            if (nod[0]== i)or(nod[1]==j):
                                Rt=Rc
                            if (nod[0]!= i)and(nod[1]!=j):
                                Rt=Rdia
                            if (cont>0):
                            	if (((Y[i,j]-Y[nod[0],nod[1]])/Rt)>(val)):
                                	val=((Y[i,j]-Y[nod[0],nod[1]])/Rt)
                                	i1=nod[0]
                                	j1=nod[1]
    val=-100
    temp=[i,j]
    """ generacion de ruta con datos del LCC====================================================="""             
    ruta=[ruta,[i,j]]
    xaxis.append(j)
    yaxis.append(i)
    temp=[i,j]
""" imprecion de coordenadas de ruta y tiempodel proceso========================================"""   
print (ruta)
tiempo1=time.clock()
tiempo = tiempo1-tiempo0
print("proceso =", tiempo ,"segundos")
"""impresion del mapa y propuesta para generar la ruta"""
plt.plot(xaxis,yaxis,'b')
plt.imshow(mapa)
plt.show()
"""FIN==========================================================================================="""