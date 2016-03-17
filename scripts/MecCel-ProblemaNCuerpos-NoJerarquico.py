
# coding: utf-8

# In[ ]:

# El problema no jerarquico de 3 cuerpos


# In[ ]:

from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
from numpy.linalg import norm
get_ipython().magic('matplotlib inline')

# Ecuaciones diferenciales para el sistema de N-cuerpos
def eom(y,t,masas):
    M=len(y);N=M/6
    r=zeros((N,3));v=zeros((N,3))
    drdt=zeros((N,3));dvdt=zeros((N,3))    
    for i in xrange(N):
        r[i]=y[3*i:3*i+3];
        v[i]=y[3*N+3*i:3*N+3*i+3]

        
    # Derivadas
    for i in xrange(N):
        drdt[i]=v[i]
        for j in xrange(N):
            if i==j:continue
            dvdt[i]+=-masas[j]/norm(r[i]-r[j])**3*(r[i]-r[j])

    # Devuelve derivadasç
    dydt=array([])
    for i in xrange(N):dydt=concatenate((dydt,drdt[i]))
    for i in xrange(N):dydt=concatenate((dydt,dvdt[i]))
    return dydt

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


# In[ ]:

# Definición del Sistema


# In[ ]:

sistema=[
    # Particula 0
    dict(
        m=10.0,
        r=[-1,0,0],
        v=[0,0,0]
    ),
    # Particula 1
    dict(
        m=5.0,
        r=[2,0,0],
        v=[0,2.2,0]
    ),
    # Particula 2
    dict(
        m=0.01,
        r=[3.6,0,0],
        v=[0,3.75133030951881,0]
    )
]
colors=['b','g','r']

# Prepara el Sistema de Partículas
Ntot=len(sistema)
masas=[]
rs=[];vs=[];ys=[]
for i in xrange(Ntot):
    particula=sistema[i]
    if particula['m']>0:
        masas+=[particula['m']]
        rs+=particula['r'];vs+=particula['v']

ys=rs+vs
M=len(ys)
N=M/6
Masa=sum(masas)


# In[ ]:

# Solución al movimiento respecto al centro de masa


# In[ ]:

Nt=1000
ts=linspace(0,20,Nt)

solucion=odeint(eom,ys,ts,args=(masas,))

rs=zeros((N,Nt,3))
vs=zeros((N,Nt,3))
for i in xrange(N):
    n=3*i
    rs[i]=solucion[:,n:n+3]
    m=3*N+3*i
    vs[i]=solucion[:,m:m+3]
    
# Posición del Centro de Masa
R = zeros((Nt,3))
for it in xrange(Nt):
    for n in xrange(N):
        R[it]+=masas[n]*rs[n,it]/Masa

# Velocidad del Centro de Masa
V=zeros((Nt,3))
for it in xrange(Nt):
    for n in xrange(N):
        V[it]+=masas[n]*vs[n,it]/Masa

# Refiere la posición de las partículas al centro de masa
for n in xrange(N):
    rs[n,:]=rs[n,:]-R
    vs[n,:]=vs[n,:]-V


# In[ ]:

# Gráfica del Movimiento


# In[ ]:

fig2d = plt.figure(figsize=(6,6))
ax2d=fig2d.gca()
for i in xrange(N):
    r=rs[i,:]
    ax2d.plot(r[:,0],r[:,1],color=colors[i])
    ax2d.plot(r[0,0],r[0,1],'o',color=colors[i],markersize=5,markeredgecolor='none')
    ax2d.plot(r[-1,0],r[-1,1],'s',color=colors[i],markersize=10,markeredgecolor='none')
    
ext=6
ax2d.set_xlim(-ext,ext)
ax2d.set_ylim(-ext,ext)
ax2d.grid()


# In[ ]:

# Problema relativo al cuerpo 2


# In[ ]:

fig2d = plt.figure(figsize=(6,6))
ax2d=fig2d.gca()
for i in xrange(N):
    r=rs[i,:]-rs[1,:]
    ax2d.plot(r[:,0],r[:,1],color=colors[i])
    ax2d.plot(r[0,0],r[0,1],'o',color=colors[i],markersize=5,markeredgecolor='none')
    ax2d.plot(r[-1,0],r[-1,1],'s',color=colors[i],markersize=10,markeredgecolor='none')
    
ext=6
ax2d.set_xlim(-ext,ext)
ax2d.set_ylim(-ext,ext)
ax2d.grid()


# In[ ]:

# Transformación al sistema rotante del cuerpo 1 y 2


# In[ ]:

r12=rs[1,:]-rs[0,:]
tetas=arctan2(r12[:,1],r12[:,0])
npos=len(tetas)

rsr=zeros_like(rs)
for i in xrange(npos):
    Rtot=rotationMatrix(tetas[i],'z')
    for n in xrange(N):
        rsr[n,i]=Rtot.dot(rs[n,i])


# In[ ]:

#  Gráfica del movimiento
fig2d = plt.figure(figsize=(6,6))
ax2d=fig2d.gca()
for i in xrange(N):
    r=rsr[i,:]
    ax2d.plot(r[:,0],r[:,1],color=colors[i])
    ax2d.plot(r[0,0],r[0,1],'o',color=colors[i],markersize=5,markeredgecolor='none')
    
ext=4
ax2d.set_xlim(-ext,ext)
ax2d.set_ylim(-ext,ext)
ax2d.grid()


# In[ ]:



