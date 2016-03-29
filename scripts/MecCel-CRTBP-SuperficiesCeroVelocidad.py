
# coding: utf-8

# # El Problema Restringido de los 3 Cuerpos: Propiedades Globales

# In[1]:

from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
from numpy.linalg import norm
get_ipython().magic(u'matplotlib inline')


# # Propiedades del Sistema

# In[3]:

mu=0.1
x1=mu
x2=1-mu
C=4.0


# # Posible Región de Movimiento

# In[4]:

x=linspace(-2,2,100)
y=linspace(-2,2,100)
X,Y=meshgrid(x,y)


# # Evaluación de superficie de Cero Velocidad

# In[5]:

Z=X**2+Y**2+2*(1-mu)/sqrt((X+x1)**2+Y**2)+2*mu/sqrt((X-x2)**2+Y**2)-C


# # Gráfico de la superficie

# In[8]:

fig=plt.figure()
ax=Axes3D(fig)

ax.plot_surface(X,Y,Z)


# # Contornos

# In[20]:

fig=plt.figure()
ax=fig.gca()

ax.plot([-x1,x2],[0,0],'o')
c=ax.contour(X,Y,Z,levels=[0,3,4,5])
ax.clabel(c)
c=ax.contourf(X,Y,Z,levels=[0,3,4,5])


# In[ ]:



