# -*- coding: utf-8 -*-
"""
Spyder Editor

Kamila Larripa
March 12 2020
Sensitivity Analysis for Macrophage Polarization
Quadstable Baseline Parameters
"""

# import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from SALib.sample import saltelli
from SALib.analyze import sobol
# to plot multiple plots
import pylab
import time
import random
import matplotlib.pyplot as plt

"""
Solve ODE for fixed parameters and initial conditions
"""

# solve the system dy/dt = f(y, t)

# model equations
def f(y, t, paras):

    x1 = y[0]
    x2 = y[1]
    

    try:
        # 16 parameter values
        a1 = paras['a1'].value
        a2 = paras['a2'].value
        b1 = paras['b1'].value
        b2 = paras['b2'].value
        k1 = paras['k1'].value
        k2 = paras['k2'].value
        l1 = paras['l1'].value
        l2 = paras['l2'].value
        n1 = paras['n1'].value
        n2 = paras['n2'].value
        p1 = paras['p1'].value
        p2 = paras['p2'].value
        q1 = paras['q1'].value 
        q2 = paras['q2'].value 
        s1 = paras['s1'].value
        s2 = paras['s2'].value
        

    except:
        a1, a2, b1, b2, k1, k2, l1, l2, n1, n2, p1, p2, q1, q2, s1, s2 = paras
    
    # the model equations 
    f0 = (a1 * (x1**n1/(x1**n1+k1**n1)) + s1) * (p1**l1 /(p1**l1 + x2**l1)) + b1 - q1 * x1
    f1 = a2 * (x2**n2/(x2**n2+k2**n2)) + s2 * (p2**l2/(p2**l2 + x1**l2)) + b2 - q2 * x2
    return [f0, f1]

# initial conditions
x10 = 2              # initial STAT1
x20 = 3                # initial STAT6

y0 = [x10, x20]     # initial condition vector
t = np.linspace(0, 10., 1000)         # time grid, 0 to 10 in 1000 steps

# baseline parameter values
a1 = 15
a2 = 8
b1 = .05
b2 = .05
k1 = 1
k2 = 1
l1 = 1
l2 = 1
n1 = 22
n2 = 6
p1 = 1
p2 = .5
q1 = 5.8
q2 = 5.8
s1 = 5 
s2 = 5



""" Global Sensitivity Analysis--- vary all parameters at the same time."""

# 16 parameters are a1, a2, b1, b2, k1, k2, l1, l2, n1, n2, p1, p2, q1, q2, s1, s2
problem = {
    'num_vars': 16, # number of parameters
    'names': ['a1', 'a2', 'b1','b2', 'k1', 'k2', 'l1', 'l2', 'n1', 'n2', 'p1', 'p2', 'q1', 'q2', 's1', 's2'],
    'bounds': [[.85*a1,1.15*a1], # vary each parameter 15% in each direction
               [.85*a2,1.15*a2],
               [.85*b1, 1.15*b1],
               [.85*b2, 1.15*b2],
               [.85*k1,1.15*k1],
               [.85*k2,1.15*k2],
               [.85*l1,1.15*l1],
               [.85*l2, 1.15*l2],
               [.85*n1,1.15*n1],
               [.85*n2,1.15*n2],
               [.85*p1,1.15*p1],
               [.85*p2, 1.15*p2],
               [.85*q1, 1.15*q1],
               [.85*q2,1.15*q2],
               [.85*s1,1.15*s1],
               [.85*s2,1.15*s2]]            
}

# generate samples
param_values = saltelli.sample(problem, 10000) # second argument is related to the number of samples to take, samples are put in rows, one column for each parameter.  10000 gives 300000 samples.
Y = np.zeros([param_values.shape[0]]) # make array to hold solution 

for i in range(param_values.shape[0]):
   
# solve the ODEs
  # set parameter values based on sample
  # 16 parameters are a1, a2, b1, b2, k1, k2, l1, l2, n1, n2, p1, p2, q1, q2, s1, s2


    a1 = param_values[i,0]     
    a2 = param_values[i,1]  
    b1 = param_values[i,2]  
    b2 = param_values[i,3]      
    k1 = param_values[i,4] 
    k2 = param_values[i,5]  
    l1 = param_values[i,6] 
    l2 = param_values[i,7] 
    n1 = param_values[i,8] 
    n2 = param_values[i,9] 
    p1 = param_values[i,10] 
    p2 = param_values[i,11] 
    q1 = param_values[i,12] 
    q2 = param_values[i,13] 
    s1 = param_values[i,14] 
    s2 = param_values[i,15] 

    soln = odeint(f, y0, t, args=((a1, a2, b1, b2, k1, k2, l1, l2, n1, n2, p1, p2, q1, q2, s1, s2), )) # solve system of ODE
    # extract each state variable, and store solution in its own matrix
    x1 = soln[:, 0]
    x2 = soln[:, 1]
    

    # let the ratio of x1 and x2 at the last time step (since we have reached steady state)
    # be the outcome of interest
    Y[i] = x1[-1]/x2[-1]
    


    
# perform analysis
# sobol.analyze will compute first, second, and total-order indices.
# Si is a Python dictionary with the keys "S1", "S2", "ST", "S1_conf", "S2_conf", and "ST_conf

Si = sobol.analyze(problem, Y)
# see results on console for 1st order and total sensitivities
print(Si['S1'])  
print(Si['ST'])

"""
plot bar graphs with error bars
"""
        
 # Si.get('S1') # returns array with values for S1    
A = Si.get('S1')
# change array to list for plotting
S1_values = A.tolist() # list of S1 values
# get ST values
B = Si.get('ST')
# change array to list for plotting
ST_values = B.tolist()
 # get confidence intervals for ST
C = Si.get('ST_conf')
# change array to list for plotting
STconf_values = C.tolist()
# get confidence intervals for S1
D = Si.get('S1_conf')
# change array to list for plotting
S1conf_values = D.tolist()
 
 # width of the bars
barWidth = 0.3
 
# Choose the height of the blue bars-- ST
bars1 = ST_values
 
# Choose the height of the cyan bars-- S1
bars2 = S1_values
 
# Choose the height of the error bars (bars1)-- CI for ST
yer1 = STconf_values
# Choose the height of the error bars (bars2)-- CI for S1
yer2 = S1conf_values
 
# The x position of bars
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]

f = plt.figure()
 
# Create blue bars
plt.bar(r1, bars1, width = barWidth, color = 'blue', edgecolor = 'black', yerr=yer1, capsize=7, label='Total Sensitivity')
 
# Create cyan bars
plt.bar(r2, bars2, width = barWidth, color = 'cyan', edgecolor = 'black', yerr=yer2, capsize=7, label='1st order Sensitivity')
 
## general layout
plt.xticks([r + barWidth for r in range(len(bars1))], ['a1', 'a2', 'b1', 'b2', 'k1', 'k2', 'l1', 'l2', 'n1', 'n2', 'p1', 'p2', 'q1', 'q2', 's1', 's2'],rotation='vertical')
plt.ylabel('Sensitivity Index')
plt.legend()
# plt.title('Sensitivity with Respect to the Ratio of Stat 1 and Stat 6 level at 10 days') # instead, add figure caption in Latex


# view and save figure (saving will allow you to put into your paper)
plt.show()
f.savefig("QuadstableSensitivityMacrophage.pdf",bbox='tight')


