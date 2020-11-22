# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:02:01 2020

@author: lenovo
"""

import math
from matplotlib import pyplot as plt  
import numpy as np  

x0=-math.pi
x1=math.pi
n=1000
dx=2*math.pi/n
h=dx*5
ad=5/(4*h)
t=[]
y=[]
w1=0
for p in range(1,1000):
    t.append(x0+dx*(p-1))
    
for i in range(1,1000):
    for j in range(i-5,i+5):
        r=(i-j)*dx
#        print(r)
        R=math.sqrt(r**2)/h 
#        print(R)
        if i==j:
            k=0
        else:
            k=r/math.sqrt(r**2)
#if(i==j) k=0.0
        if R>1:
            w=0
        else:
            w=ad*(1-R)**2*(-12*R)*(1/h)*k
            print(w) 
        w1=math.sin(x0+dx*(j-1))*w*dx+w1
        print(w1)
    
    y.append(w1)
    w1=0
    
plt.scatter(t,y,s=10)
plt.show