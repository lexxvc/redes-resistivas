import numpy as np                # funciones numéricas (arrays, matrices, etc.)
import matplotlib.pyplot as plt   # funciones para representación gráfica
img='face.jpg'              #integracion de imagen
I = plt.imread(img)         # se pasa la imagen a valores matriciales
plt.imshow(I,vmin=0,vmax=1) 
x=I.shape[0]                #dimensiones de la  imagen
y=I.shape[1]
size=(x,y)
data=np.ones(size)         # matris de tamaño MxN de ceros
#print (data)
#plt.show()
for i in range (x):         # primer face de llenado 
    data[x-1,i]=i+1
    data[x-i-1,0]=i+1
#print (data)
t2=np.zeros((x-1,y-1))
a=t2.size
b=[]
for k in range (a,0,-1):
    b.append(k)
h=0
for i in range (x-2,-1,-1):
        for j in range (0,y-1):
            t2[i,j]= b[h]
            h=h+1
data2=t2.transpose()

for i in range (x-1):
        for j in range (1,y):
            data[i,j]=data2[i,j-1]
print(data)