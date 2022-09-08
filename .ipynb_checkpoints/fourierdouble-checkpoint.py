#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import matplotlib.pyplot as plt
import sympy as smp
from mpl_toolkits.mplot3d import axes3d, Axes3D
from sympy.plotting import plot3d


# In[2]:


x,y,a,b,p,wm,v,d,E,f,h=smp.symbols('x y a b p wm v d E f h')
m,n,i,j=smp.symbols('m n i j', int=True)


# In[3]:




# In[4]:


def sfourier_doublesin(f,a,b,nn,mm):
    freq=smp.sin(x*smp.pi*m/a)*smp.sin(y*smp.pi*n/b)
    fmn=4/a/b*smp.integrate(smp.integrate(f*freq,(x,0,a)),(y,0,b))
    F=smp.summation((smp.summation(fmn*freq,(n,1,nn))),(m,1,mm))
    plot3d(F, (x, 0,a), (y, 0,b))
    plot3d(f, (x, 0,a), (y, 0,b))
    V1=smp.integrate(smp.integrate(F,(x,0,a)),(y,0,b))
    V2=smp.integrate(smp.integrate(f,(x,0,a)),(y,0,b))
    erro=V1/V2
    print(erro.evalf())
    return fmn,F




def sfourier_value(f,a,b,nn,mm):
    freq=smp.sin(x*smp.pi*m/a)*smp.sin(y*smp.pi*n/b)
    fmn=4/a/b*smp.integrate(smp.integrate(f*freq,(x,0,a)),(y,0,b))
    F=smp.summation((smp.summation(fmn*freq,(n,1,nn))),(m,1,mm))
    return fmn,F


def sfourier_etaxsi(f,a,b,xsi,eta,c,d,nn,mm):
    freq=smp.sin(x*smp.pi*m/a)*smp.sin(y*smp.pi*n/b)
    fmn=4/a/b*smp.integrate(smp.integrate(f*freq,(x,xsi-c/2,xsi+c/2)),(y,eta-d/2,eta+d/2))
    F=smp.summation((smp.summation(fmn*freq,(n,1,nn))),(m,1,mm))
    plot3d(F, (x, 0,a), (y, 0,b))

   


    return fmn,F