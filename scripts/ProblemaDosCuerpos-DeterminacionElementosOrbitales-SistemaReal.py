
# coding: utf-8

# # Determinación de los Elementos Orbitales de un Sistema Real

# In[78]:

from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import bisect
from scipy.special import jn
from time import time as timeit
get_ipython().magic(u'matplotlib inline')
norm=linalg.norm

DEG=pi/180
RAD=180/pi
def show(r):"["+",".join(r.tolist())+"]"


# # Definición del Sistema

# In[79]:

# Tiempo 2016-06-19 02:00:00
# Dia Juliano
t=2457558.58333

# Se utilizo la interface web del Sistema Horizons del JPL : 
# http://ssd.jpl.nasa.gov/horizons.cgi

# Tierra
m1 = 5.96e24
r1 = array([-3.164891183028434E-02,-1.012904175107438E+00,-1.277405128356421E-04])
v1 = array([1.690601584062892E-02,-6.545632536119304E-04,-6.217673594760544E-08])

# Luna
m2 = 7.34e22
r2 = array([-3.246216458663009E-02,-1.015435304053814E+00,1.050538209588863E-04])
v2 = array([1.745520310088857E-02,-8.109331213061256E-04,1.516106948850730E-06])


# # Unidades canónicas del problema

# In[80]:

# Cantidades de referencia 
AU=1.5e8*1e3 # Unidad Astronómica, m
G=6.67e-11 # m^3/(kg s)
Rt = 6.4e3*1e3 # Radio de la Tierra, m

# Unidades
ul=Rt
um=m1
ut=sqrt(ul**3/(G*um))

print "ul, ut, um = ", ul, ut, um


# # Propiedades del sistema en unidades canónicas

# In[81]:


t=t*86400/ut

m1=m1/um
r1=r1*AU/ul
v1=v1*AU/(86400)*(ut/ul)

m2=m2/um
r2=r2*AU/ul
v2=v2*AU/(86400)*(ut/ul)

print "Tierra: m, r, v: ",m1,r1,v1
print "Luna: m, r, v: ",m2,r2,v2


# # Propiedades dinámicas

# In[82]:

# Posición inicial del Centro de masa
mu = m1+m2
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

print "mu = ",mu

print "h = ",h
print "eps = ",eps
print "evec = ",evec


# # Propiedades de la cónica

# In[83]:

e=norm(evec)
p=norm(h)**2/mu

print "Orbital geometrical properties:"
print "p,e = ",p,e

if e<1:
    a=p/(1-e**2)
    n=sqrt(mu/a**3)
    P=2*pi/n
    print "a = ",a*ul
    print "Movimiento orbital medio = ",n
    print "Período = ",P*ut/86400
else:
    print "La órbita es una hipérbola"


# # Orientación en el espacio

# In[84]:

# Vectores unitarios
ax = array([1,0,0])
ay = array([0,1,0])
az = array([0,0,1])

# Vectores plano orbital
I = az
N = cross(az,h)
print N

# Elementos orbitales clásicos
i = arccos(h[2]/norm(h))
o = arccos(N[0]/norm(N))
w = arccos(dot(N,evec)/(norm(N)*norm(evec)))

# El dificil
absf = arccos(dot(evec,r)/(norm(evec)*norm(r)))
if dot(cross(r,evec),h)>0:f=-absf
else:f=absf

print "Classical orbital elements i, Omega, omega, f:",i*RAD,o*RAD,w*RAD,f*RAD


# Tiempo de paso por el periapsis
E = 2 * arctan(sqrt((1-e)/(1+e))*tan(f/2))
M = E - e*sin(E)
tp = t-M/n

print "Tiempo de paso por el periapsis (jd):",(tp*ut/86400)


# # Resumen

# In[87]:

print "# Sistema Tierra-Luna en una fecha específica:"
print "t = ",t

print "# Tierra"
print "m1 = ",m1
print "r1 = array(",r1.tolist(),")"
print "v1 = array(",v1.tolist(),")"

print "# Luna"
print "m2 = ",m2
print "r2 = array(",r2.tolist(),")"
print "v2 = array(",v2.tolist(),")"

print "# Centro de Masa"
print "rcm = array(",rcm.tolist(),")"
print "vcm = array(",vcm.tolist(),")"

print "# Unidades del sistema a SI"
print "ul = ",ul
print "ut = ",ut
print "um = ",um

print "# Propiedades dinámicas (unidades canónicas)"
print "h = array(",h.tolist(),")"
print "evec = array(",evec.tolist(),")"
print "n = ",n
print "P = ",P

print "# Elementos orbitales (en unidades canónicas y radianes)"
print "a = ",a
print "e = ",e
print "i = ",i
print "o = ",o
print "w = ",w
print "f = ",f
print "tp = ",tp


# In[ ]:



