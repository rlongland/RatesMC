#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 10:17:12 2022

@author: larswalker
"""

import matplotlib.pyplot as plt

with open('RatesMC.Lars') as f:

    lines = f.readlines()
    
    



fast = []
slow = []


for value in lines:
    dummy = value.split(' ')
    fast.append(float(dummy[0]))
    slow.append(float(dummy[1]))
    
  
difference = []
total = []
counter = 0
    
for value in fast:
    diff = value - slow[counter]
    difference.append(diff)
    total.append(counter)
    counter += 1



#Get the points from 55,000 to 60,000
def zoom(Y_Points, start, stop):
    zoomed_y = []
    zoomed_x = []
    for num in range(start, stop):
        zoomed_y.append(Y_Points[num])
        zoomed_x.append(num-start)
    
    plt.title("Zoomed")    
    plt.scatter(zoomed_x, zoomed_y)
    plt.show()

   
    
   
plt.title("Fast")    
plt.scatter(total, fast)
plt.show()


plt.title("Slow")
plt.scatter(total, slow)
plt.show()


plt.title("Fast and Slow")
plt.scatter(total, slow)
plt.scatter(total, fast)
plt.show()


zoom(slow, 0, 1000)

    
plt.title("difference")
plt.scatter(total, difference)