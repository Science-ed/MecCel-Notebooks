
# coding: utf-8

# # Determinación de los Elementos Orbitales

# In[104]:

from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import bisect
from scipy.special import jn
from time import time as timeit
get_ipython().magic(u'matplotlib inline')
norm=linalg.norm

DEG=pi/180
RAD=180/pi


# # Definición del Sistema

# In[106]:

t=0

m1 = 5.0
r1 = array([1,0,0])
v1 = array([0,1,1])

m2 = 3.0
r2 = array([-1,0,0])
v2 = array([1,-1,0])


# # Propiedades dinámicas

# In[108]:

# Posición inicial del Centro de masa
rcm = (m1*r1+m2*r2)/mu
vcm = (m1*v1+m2*v2)/mu

# Posición relativa
mu=m1+m2
r=r1-r2
v=v1-v2

# Constantes de movimiento
h=cross(r,v)
eps=0.5*norm(v)**2-mu/norm(r)
evec=-cross(h,v)/mu-r/norm(r)

print "h = ",h
print "eps = ",eps
print "evec = ",evec


# # Propiedades de la cónica

# In[110]:

e=norm(evec)
p=norm(h)**2/mu

print "Orbital geometrical properties:"
print "p,e = ",p,e

if e<1:
    a=p/(1-e**2)
    n=sqrt(mu/a**3)
    P=2*pi/n
    print "Movimiento orbital medio = ",n
    print "Período = ",P
else:
    print "La órbita es una hipérbola"


# # Orientación en el espacio

# In[115]:

# Vectores unitarios
ax = array([1,0,0])
ay = array([0,1,0])
az = array([0,0,1])

# Vectores plano orbital
I = az
N = cross(az,h)

# Elementos orbitales clásicos
i = arccos(h[2]/norm(h))
o = arccos(N[0]/norm(N))
w = arccos(dot(N,evec)/(norm(N)*norm(evec)))
f = arccos(dot(evec,r)/(norm(evec)*norm(r)))

print "Classical orbital elements i, Omega, omega, f:",i*RAD,o*RAD,w*RAD,f*RAD

# Tiempo de paso por el periapsis
E = 2 * arctan(sqrt((1-e)/(1+e))*tan(f/2))
M = E - e*sin(E)
tp = t-M/n

print "Tiempo de paso por el periapsis:",tp


# In[ ]:



