
# coding: utf-8

# # Gráfica de una cónica usando su ecuación algebraica

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')


# # Variables

# In[44]:

# Coeficientes que definen la cónica
# A x^2 + B x y + C y^2 + D x + E y + F = 0
A=1
B=0
C=0
D=0
E=1
F=0

# Matriz de Coeficientes
Cmat = [[A,B/2,D/2],[B/2,C,E/2],[D/2,E/2,F]]

# "Signo" eta de la matriz
det = np.linalg.det(Cmat)
if det<=0:eta=+1
else:eta=-1


# In[40]:

# Discriminante
d = B**2 - 4*A*C
print "Discriminante: ",d
# Excentricidad
e = np.sqrt(2*np.sqrt((A-C)**2+B**2)/(eta*(A+C)+np.sqrt((A-C)**2+B**2)))
print "Excentricidad: ",e


# # Cálculo de puntos

# In[41]:

# Vectores con puntos
xs = []
ys = []


# In[42]:

# Para un conjunto de puntos en x
for x in np.linspace(-10,10,1000):
    coeficientes=[C,E+B*x,A*x**2+D*x+F]
    yes=np.roots(coeficientes)
    for y in yes:
        if np.abs(y.imag)<1e-5:
            xs+=[x]
            ys+=[y]
            
ext=max(max(np.abs(xs)),max(np.abs(ys)))


# # Gráfica

# In[43]:

plt.figure(figsize=(6,6))
ax=plt.gca()
ax.plot(xs,ys,'o',markersize=1)
ax.set_xlim((-ext,ext))
ax.set_ylim((-ext,ext))


# In[ ]:



