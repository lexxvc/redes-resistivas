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
print ('a',V)
for i in range (x):
    V.append([i,0])
    
print ('b',V)
for i in range (y-1):
    V.append([0,i+1])
print ('c',V)

for i in range(1,x):
    for j in range (1,y):
        V.append([i,j])
print ('d',V)
nopsV = V.__len__()
print(nopsV)

"""eval tipo 0"""
L=[1]
"""eval tipo 1"""
k = 0
for i in range (y-1):
    L.append(k)
    L[k+1]=(1/4)*abs(abs(x2-i)+x2-i)-(1/4)*abs(x2-i-abs(x2-i))+(1/4)*abs(x2-i)+(1/4)*abs(x2-i)
    k=k+1
for i in range (x-1,y-1):
    L.append(k)
    L[k+1]=(1/4)*abs(abs(x2-i)+x2-i)-(1/4)*abs(x2-i-abs(x2-i))+(1/4)*abs(x2-i)+(1/4)*abs(x2-i)
    k=k+1
print(L)
"""eval tipo 2"""
for i in range (x-1):
    for j in range (y-1):
        L.append(k)
        L[k+1]=(1/4)*abs(abs(x2-i)+x2-i)-(1/4)*abs(x2-i-abs(x2-i))+(1/4)*abs(x2-i)+(1/4)*abs(x2-i)
        k=k+1
print(L)
nopsL = L.__len__()
print(nopsL)
"""###############################"""
lamda=[]
for i in range (x*y): #deberia ser x*y pero se sale de los limites
    lamda.append(L[i])

A=[]
for i in range (x*y):
    A=A ,eval(lamda,[x1=V[i,0],x2=V[i,2]])
print(A) 
    