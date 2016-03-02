
# coding: utf-8

# In[48]:

from numpy import *


# In[49]:

#Calcula la fuerza gravitacional de la Tierra y el Sol sobre la Luna (unidades SI)

Msol = 2e30 # kg
Mtierra = 6e24 # kg
Mluna = 7.34e22 # kg
AU = 1.5e11 # m

Rts = AU # m
Rlt = 3.8e8 # m

Gconst = 6.67e-11 # m^3 / kg s^2 

Rls = Rts + Rlt

als = Gconst * Msol / Rls**2
alt = Gconst * Mtierra / Rlt**2
a = als+alt

print "Aceleración hacia el Sol: ",als
print "Aceleración hacia la Tierra: ",alt
print "Aceleración Total: ",a


# In[50]:

#Calcula lo mismo unidades canónicas
G = 1 # uL^3 / uM uT^2
uM = Msol
uL = AU
uT = sqrt(uL**3/(Gconst*uM))

print "Unidades: uL = %.2e m, uM = %.2e kg, uT = %.2e s"%(uL,uM,uT)


# In[51]:

Msol = 1.0
Mtierra = Mtierra / uM
Rts = 1.0
Rlt = Rlt/uL
Rls = Rts+Rlt

als = Msol/Rls**2
alt = Mtierra/Rlt**2
a = als+alt

print "Aceleración hacia el Sol: ",als
print "Aceleración hacia la Tierra: ",alt
print "Aceleración Total: ",a


# In[ ]:



