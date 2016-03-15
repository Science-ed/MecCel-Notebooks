
# coding: utf-8

# # "Brujería" con el problema de los dos cuerpos

# In[172]:

from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import bisect
from scipy.special import jn
from time import time as timeit
get_ipython().magic(u'matplotlib inline')
norm=linalg.norm

DIA=86400.0 # s
DEG=pi/180
RAD=180/pi
AU=1.5e11 # m

# Solución ultraeficiente y precisa a la ecuación de Kepler
def Esol(M,e):
    """
    Mikkola, 1991
    Code at: http://smallsats.org/2013/04/20/keplers-equation-iterative-and-non-iterative-solver-comparison/
    """
    Ecorr=0;Esgn=1.0
    if M>180*DEG:
        M=360.0*DEG-M
        Ecorr=360*DEG
        Esgn=-1.0
    if e==0:return Ecorr+Esgn*M
        
    a=(1-e)*3/(4*e+0.5);
    b=-M/(4*e+0.5);
    y=(b**2/4 +a**3/27)**0.5;
    x=(-0.5*b+y)**(1./3)-(0.5*b+y)**(1./3);
    w=x-0.078*x**5/(1 + e);
    E=M+e*(3*w-4*w**3);

    #NEWTON CORRECTION 1
    sE=sin(E)
    cE=cos(E)

    f=(E-e*sE-M);
    fd=1-e*cE;
    f2d=e*sE;
    f3d=-e*cE;
    f4d=e*sE;
    E=E-f/fd*(1+                  f*f2d/(2*fd*fd)+                  f*f*(3*f2d*f2d-fd*f3d)/(6*fd**4)+                  (10*fd*f2d*f3d-15*f2d**3-fd**2*f4d)*                  f**3/(24*fd**6))

    #NEWTON CORRECTION 2
    f=(E-e*sE-M);
    fd=1-e*cE;
    f2d=e*sE;
    f3d=-e*cE;
    f4d=e*sE;
    E=E-f/fd*(1+                  f*f2d/(2*fd*fd)+                  f*f*(3*f2d*f2d-fd*f3d)/(6*fd**4)+                  (10*fd*f2d*f3d-15*f2d**3-fd**2*f4d)*                  f**3/(24*fd**6))
    
    E=Ecorr+Esgn*E
    return E

# Convierte una posiciones devueltas por Horizons en vectores de posición y velocidad
def str2vec(cadena):
    lista=cadena.split("\n")
    pr=lista[1].split(" ")
    pv=lista[2].split(" ")
    r=[]
    for x in pr:
        if x=='':continue
        r+=[float(x)]
    r=array(r)
    v=[]
    for x in pv:
        if x=='':continue
        v+=[float(x)]
    v=array(v)
    return r,v


# # Definición del Sistema

# In[173]:

# Esto viene de MecCel-ProblemaDosCuerpos-DeterminacionElementosOrbitales-SistemaReal
# Sistema Tierra-Luna en una fecha específica:
t =  263618127.507
# Tierra
m1 =  1.0
r1 = array( [-743.1501836105128, -23784.06966240977, -2.9994833970080177] )
v1 = array( [3.700732439155043, -0.14328411193718754, -1.3610507990244552e-05] )
# Luna
m2 =  0.0122932772082
r2 = array( [-762.2462251566171, -23843.50326794187, 2.4667756905271423] )
v2 = array( [3.8209497114191127, -0.17751352750957722, 0.00033187631076013165] )
# Centro de Masa
rcm = array( [-743.3820857060726, -23784.79142339025, -2.933101213338583] )
vcm = array( [3.70219235624274, -0.14369979354229206, -9.414920278357005e-06] )
# Unidades del sistema a SI
ul =  6371000.0
ut =  805.456982823
um =  5.97237e+24
# Propiedades dinámicas (unidades canónicas)
mu =  1.01229327721
h = array( [0.16657332663144261, 0.6637361876366484, 7.798592279757015] )
evec = array( [-0.0408065065347476, -0.02235168365874296, 0.0027739489467564493] )
n =  0.00212887649751
P =  2951.40902468
# Elementos orbitales (en unidades canónicas y radianes)
a =  60.6738852571
e =  0.0466096935273
i =  0.0875248484444
o =  2.89570769884
w =  0.748909787804
f =  -2.38445516812
tp =  263619216.709


# # Fecha en la que se quiere calcular la posición

# In[187]:

# 2016-07-01 02:00:00 en dias Julianos
tfuturo = 2457570.583333330 * 86400 / ut

# 2016-12-01 02:00:00 en dias Julianos
# tfuturo = 2457723.583333330 * 86400 / ut


# # Cálculo de la posición en el plano fundamental

# In[188]:

# Anomalías
M=mod(n*(tfuturo-tp),2*pi)
E=Esol(M,e)
f=2*arctan(sqrt((1+e)/(1-e))*tan(E/2))

print "Anomalía verdadera: ",f*RAD

# r,v
r=a*(1-e*cos(E))
v=sqrt(mu*(2/r-1/a))
print "r,v = ",r,v

# Angulo entre la velocidad y el vector posición
alpha=arcsin(norm(h)/(norm(r)*norm(v)))
print "Angulo velocidad radio vector: ",alpha*RAD

# Ángulo de la velocidad
teta = f+alpha
print "Angulo de la velocidad: ",teta*RAD

# Coordenadas en el plano fundamental
x = r*cos(f)
y = r*sin(f)
z = 0
rpp = array([x,y,z])

vx = v*cos(teta)
vy = v*sin(teta)
vz = 0
vpp = array([vx,vy,vz])

print "Posición en el plano fundamental: ",rpp.tolist()
print "Velocidad en el plano fundamental: ",rpp.tolist()


# # Posición de la Tierra y la Luna en el plano fundamental

# In[189]:

r1pp=+m2/mu*rpp
v1pp=+m2/mu*vpp

r2pp=-m1/mu*rpp
v2pp=-m1/mu*vpp


# # Matriz de rotación al sistema de referencia de observación

# In[190]:

# Matriz de Rotación
def rotationMatrix(t,axis):
    R=identity(3)
    r=array([[cos(t),sin(t)],[-sin(t),cos(t)]])
    if axis=='z':R[0:2,0:2]=r
    elif axis=='x':R[1:3,1:3]=r
    else:
        R[0,0]=r[0,0];R[0,2]=r[0,1]
        R[2,0]=r[1,0];R[2,2]=r[1,1]
    return R

Rtot=(rotationMatrix(-o,'z').dot(rotationMatrix(-i,'x'))).dot(rotationMatrix(-w,'z'))


# # Posición y velocidad en el sistema de referencia de la observación

# In[191]:

r1 = Rtot.dot(r1pp)
v1 = Rtot.dot(v1pp)

r2 = Rtot.dot(r2pp)
v2 = Rtot.dot(v2pp)

print "Posición de la Tierra r, v:",r1,v1
print "Posición de la Luna r, v:",r2,v2


# # En unidades del sistema Horizons del JPL

# In[192]:

r1=r1*ul/AU
v1=v1*ul/ut*(DIA/AU)

r2=r2*ul/AU
v2=v2*ul/ut*(DIA/AU)

print "Posición y velocidad de la Tierra:"
print "%+.7e %+.7e %+.7e "%(r1[0],r1[1],r1[2])
print "%+.7e %+.7e %+.7e]"%(v1[0],v1[1],v1[2])

print "Posición y velocidad de la Luna:"
print "%+.7e %+.7e %+.7e "%(r2[0],r2[1],r2[2])
print "%+.7e %+.7e %+.7e]"%(v2[0],v2[1],v2[2])


# # Valores calculados con Horizons

# In[196]:

# 2016-07-01 02:00:00
luna="""
1.420963124590351E-03  1.944834162289680E-03 -1.991995015794231E-04
-4.959947218313727E-04  3.586965896461778E-04 -1.919821623749072E-05
"""
  
tierra="""
-1.747789887292984E-05 -2.392153197002750E-05  2.450161220859895E-06
6.100753383149802E-06 -4.411981290301303E-06  2.361387682289336E-07
"""

r1n,v1n=str2vec(luna)
r2n,v2n=str2vec(luna)


# # Errores

# In[197]:

print "Errores en posición de la Luna (%):",(r2-r2n)/norm(r2n)*100
print "Errores en velocidad de la Luna (%):",(v2-v2n)/norm(v2n)*100


# # Igual pero en la fecha de diciembre

# In[198]:

# 2016-12-01 02:00:00
luna="""
-1.711627476952467E-04 -2.647728759600073E-03  2.193944996516806E-04
5.578391890816546E-04 -2.411335140327116E-05 -1.672378277158184E-05
"""

tierra="""
2.105308113391589E-06  3.256716145718358E-05 -2.698560442442569E-06
-6.861442612690692E-06  2.965951121604893E-07  2.057031452795932E-07
"""
r1n,v1n=str2vec(luna)
r2n,v2n=str2vec(luna)

print "Errores en posición de la Luna (%):",(r2-r2n)/norm(r2n)*100
print "Errores en velocidad de la Luna (%):",(v2-v2n)/norm(v2n)*100


# In[ ]:



