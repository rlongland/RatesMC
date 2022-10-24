#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 10:17:12 2022

@author: Lars Walker
email: larsw00@outlook.com
"""

import matplotlib.pyplot as plt
import numpy


with open('RatesMC.Lars') as f:

    lines = f.readlines()

#Initiallize all the lists for each temp (f for fast, s for slow):
#Initialize List of List of Lists
LLL = []   
Temps = []
TL = [] 
R =[]
g0 = []
g1 = []
g2 = []




for value in lines:
    dummy = value.split(' ')
    T = float(dummy[2])
    if T in Temps:
        R.append(float(dummy[0])/float(dummy[1]))
        g0.append(float(dummy[3]))
        g1.append(float(dummy[4]))
        g2.append(float(dummy[5]))
    else: 
        Temps.append(T)
        if (bool(R)):
            TL.extend([R,g0,g1,g2])
            LLL.insert(Temps.index(T), TL)
            TL = [] 
            R =[]
            g0 =[]
            g1 = []
            g2 = []
            
            
TL.extend([R,g0,g1,g2])
LLL.insert(Temps.index(T), TL) 
print("Possible Tempuratures: ")
for a in Temps:
    print(a)

        



def total(f):
    tot =[]
    for num in range(0,len(f)):
        tot.append(num)
    return tot



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

    
def Plot(T, index):
    # num of resonances (constant for this input file), res energy, gamma 1, gamma 2, gamma 3
    X = []
    Y = []
    title =f'at T = {T}'
    Y = LLL[Temps.index(T)][0]
    
    
        
    if (index == 0):
        plt.title("No Variable "+title)
        X = [total(Y)]
        
    else:
        plt.title(f' G{index-1} '+title)
        X = LLL[Temps.index(T)][index]
        
        
    plt.scatter(X,Y)
    plt.show()


Plot(.001,1)
Plot(0.001,0)


