
# coding: utf-8

# In[146]:

import numpy as np
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')


# In[147]:

# Dibujando una elipse a partir de su ecuación en coordenadas polares
# r = p / (1 + e cos f)


# In[148]:

p = 1
e = 0.2


# In[149]:

f=np.linspace(-np.pi,np.pi,1000)
rs=p/(1+e*np.cos(f))
xs=rs*np.cos(f)
ys=rs*np.sin(f)


# In[150]:

fig=plt.figure(figsize=(6,6))
ax=fig.gca()

ax.plot(xs,ys,'-')
ax.grid()
ext=1.5
ax.set_xlim((-ext,ext))
ax.set_ylim((-ext,ext))


# In[151]:

# ¿Qué pasa en el caso de una hipérbola? Debemos usar el máximo valor de la anomalía verdadera
e=1.5
fmax=np.arccos(-1/e)

f=np.linspace(-0.8*fmax,0.8*fmax,1000)
rs=p/(1+e*np.cos(f))
xs=rs*np.cos(f)
ys=rs*np.sin(f)

fig=plt.figure(figsize=(6,6))
ax=fig.gca()

ax.plot(xs,ys,'-')
ax.grid()
ext=0.8
ax.set_xlim((-ext,ext))
ax.set_ylim((-ext,ext))


# In[152]:

# Dibujando una cónica rotada
p=1
e=0.8
theta=100.0*np.pi/180

fs=np.linspace(-np.pi,np.pi,1000)
rs=p/(1+e*np.cos(fs))
xs=rs*np.cos(fs)
ys=rs*np.sin(fs)


# In[153]:

# La rotación se logra aplicando una matriz de rotación
sint=np.sin(theta)
cost=np.cos(theta)

Mrot = np.array([[cost,-sint],[sint,cost]])


# In[154]:

# Las coordenadas rotadas son
xps=[]
yps=[]
for x,y in zip(xs,ys):
    r=[x,y]
    rp=Mrot.dot(r)
    xps+=[rp[0]]
    yps+=[rp[1]]


# In[155]:

fig=plt.figure(figsize=(6,6))
ax=fig.gca()

ax.plot(xps,yps,'-')
ax.grid()
ext=5
ax.set_xlim((-ext,ext))
ax.set_ylim((-ext,ext))


# In[156]:

# Anomalía Excéntrica


# In[163]:

# Podemos repetir gráficos similares usando la anomalía excéntrica
a=1
b=0.5
Es=np.linspace(0,2*np.pi,100)

xs=a*np.cos(Es)
ys=b*np.sin(Es)


# In[164]:

fig=plt.figure(figsize=(6,6))
ax=fig.gca()

ax.plot(xs,ys,'-')
ax.grid()

ext=1.5
ax.set_xlim((-ext,ext))
ax.set_ylim((-ext,ext))


# In[159]:

# En coordenadas polares
a=1
e=0.5
Es=np.linspace(0,2*np.pi,100)
fs=2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(Es/2))

rs=a*(1-e*np.cos(Es))
xs=rs*np.cos(fs)
ys=rs*np.sin(fs)


# In[165]:

fig=plt.figure(figsize=(6,6))
ax=fig.gca()

ax.plot(xs,ys,'-')
ax.grid()

ext=1.5
ax.set_xlim((-ext,ext))
ax.set_ylim((-ext,ext))


# In[166]:

# O bien usando las propiedades de la anomalía excéntrica
b=a*np.sqrt(1-e**2)
xs=a*np.cos(Es)-a*e
ys=b*np.sin(Es)


# In[167]:

fig=plt.figure(figsize=(6,6))
ax=fig.gca()

ax.plot(xs,ys,'-')
ax.grid()

ext=1.5
ax.set_xlim((-ext,ext))
ax.set_ylim((-ext,ext))


# In[206]:

# La hipérbola también se puede graficar usando la anomalía excéntrica
a=1
e=1.5
fmax=0.999*np.arccos(-1/e)
Fmax=2*np.arctanh(np.sqrt((e-1)/(e+1))*np.tan(fmax/2))

Fs=np.linspace(-Fmax,Fmax,1000)

sinhF=np.sinh(Fs)
rs=a*np.sqrt(1+e**2*sinhF**2)

xs=rs*np.cosh(Fs)
ys=rs*np.sinh(Fs)


# In[207]:

fig=plt.figure(figsize=(6,6))
ax=fig.gca()

ax.plot(xs,ys,'-')
ax.grid()

ext=3.5
ax.set_xlim((-ext,ext))
ax.set_ylim((-ext,ext))


# In[ ]:



