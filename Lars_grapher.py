#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 10:17:12 2022

@author: larswalker
"""

import matplotlib.pyplot as plt
import numpy

with open('RatesMC.Lars') as f:

    lines = f.readlines()

#Initiallize all the lists for each temp (f for fast, s for slow):

f1 =[]
f05 = [] 
f01 = []
f001 = []
s1 =[]
s05 =[]
s01 = []
s001 = []




for value in lines:
    dummy = value.split(' ')
    T = float(dummy[2])
    if ( T == .1):
        f1.append(float(dummy[0]))
        s1.append(float(dummy[1]))
        
    elif ( T == .05):
        f05.append(float(dummy[0]))
        s05.append(float(dummy[1]))
        
    elif ( T == .01):
        f01.append(float(dummy[0]))
        s01.append(float(dummy[1]))
        
    elif ( T == .001):
        f001.append(float(dummy[0]))
        s001.append(float(dummy[1]))
        
    
    
    
  
def total(f):
    T = []
    for num in range(0,len(f)):
        T.append(num)
        
    print (len(T))
    return T



#Get the points from 55,000 to 60,000
def zoom(Y_Points, start, stop):
    zoomed_y = []
    zoomed_x = []
    for num in range(start, stop):
        zoomed_y.append(Y_Points[num])
        zoomed_x.append(num-start)
    
     
    plt.title(f"Zoomed from {start} to {stop}")
        
        
    plt.scatter(zoomed_x, zoomed_y)
    plt.show()
    

def get_ratio(f,s):
    R = []
    
    for num in range(0,len(f)):
        ratio = f[num]/s[num]
        R.append(ratio)
    return R

    
def Plot(T, index):
    # num of resonances (constant for this input file), res energy, gamma 1, gamma 2, gamma 3
    X = []
    Y = []
    title =""
    
    
    if ( T == .1):
        Y= get_ratio(f1, s1)
        title = "Ratio at T  = .1"
        
    elif ( T == .05):
        Y = get_ratio(f05, s05)
        title = "Ratio at T  = .05"
        
    elif ( T == .01):
        Y = get_ratio(f01, s01)
        title = "Ratio at T  = .01"
        
    elif ( T == .001):
        Y = get_ratio(f001, s1)
        title = "Ratio at T  = .001"
        
    if (index == 0):
        plt.title("No Variable "+title)
        X = total(Y)
        
        
        
    plt.scatter(X,Y)
    plt.show()
    




Plot(.1,0)
Plot(.05,0)
Plot(.01,0)
Plot(.001,0)
