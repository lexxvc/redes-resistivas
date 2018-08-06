import numpy as np               
import matplotlib.pyplot as plt 
from sympy import *
x1=Symbol('x1')
x2=Symbol('x2')
img='face.jpg'              
I = plt.imread(img)       
plt.imshow(I,vmin=0,vmax=1) 
x=I.shape[0]                
y=I.shape[1]
size=(x,y)
print(size)
V=[] 
#print ('a',V)
for i in range (x):
    V.append([i,0])
    
#print ('b',V)
for i in range (y-1):
    V.append([0,i+1])
#print ('c',V)

for i in range(1,x):
    for j in range (1,y):
        V.append([i,j])
print ('Datos ',V)
nopsV = V.__len__()
print(nopsV)

"""eval tipo 0"""
L=[1]

"""eval tipo 1"""
print('tipo 1')
k = 0
l = 0
#print('L',k,' =',L)
for i in range (x-1+y-1):
    if (i < x-1):
        L.append(k)
        L[k+1]=(1/4)*abs(abs(x1-k)+x1-k)-(1/4)*abs(x1-k-abs(x1-k))+(1/2)*abs(x1-k)
        k=k+1
        #print('AL',k,' =',L[k])
    if (i>=x-1):
        L.append(k)
        L[k+1]=(1/4)*abs(abs(x2-l)+x2-l)-(1/4)*abs(x2-l-abs(x2-l))+(1/2)*abs(x2-l)
        k=k+1
        l=l+1
        #print('BL',k,' =',L[k])    
"""eval tipo 2"""
print('tipo 2')
for i in range (x-1):
    for j in range (y-1):
        L.append(k)
        L[k+1]=(1/4)*abs(abs(x1-i)+x2-j)-(1/4)*abs(x1-i-abs(x2-j))+(1/4)*abs(x1-i)+(1/4)*abs(x2-j)-(1/4)*abs(x1-i-(x2-j))
        k=k+1
        #print(L[k])
"""###############################"""
lamda=[]
for i in range (x*y):
    lamda.append(i)
    lamda[i]=L[i]
    #print(lamda)
"""#################################"""
    


    