### CODE FOR EX 1.C
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import least_squares
import seaborn as sns;
import matplotlib.pyplot as plt

# Define parameters
theta = 0.67
delta = 1/16
beta = 4/(1-theta+4*(1-delta))
h = 0.31
z = 1
S = 56
c0 = 0.6286

# Steady state
k = (((1-theta)/( (beta)**(-1) - 1 + delta))**(1/theta))*z*h
c = (k**(1-theta))*((z*h)**theta) - delta*k
y  = (k**(1-theta))*((z*h)**theta)

# New steady state
kk = (((1-theta)/( (beta)**(-1) - 1 + delta))**(1/theta))*2*z*h
cc = (kk**(1-theta))*((2*z*h)**theta) - delta*kk
yy  = (kk**(1-theta))*((2*z*h)**theta)

# Transition
k_trans = np.zeros((1,S+1))
c_trans = np.zeros((1,S+1))
s_trans = np.zeros((1,S+1))
y_trans = np.zeros((1,S+1))
h_trans = np.zeros((1,S+1)) + 0.31 # Labor is constant
deviation = np.zeros((1,S+1))

# Initial values for transition paths
k_trans[0,0] = k
c_trans[0,0] = c0
y_trans[0,0] = y
s_trans[0,0] = y - c

# Index for loop
indexktrans = np.arange(0, S, 1).tolist()
   
# Compute transition paths
for i in indexktrans:
    k_trans[0,i+1] = (k_trans[0,i]**(1-theta))*((2*z*h)**theta) + (1-delta)*k_trans[0,i] - c_trans[0,i] # from law of motion
    c_trans[0,i+1] = beta*c_trans[0,i]*(1-delta +(1-theta)*(k_trans[0,i+1]**(-theta))*((2*z*h)**theta) )   #from euler equation 
    y_trans[0,i+1] = (k_trans[0,i+1]**(1-theta))*((2*z*h)**theta)   #from production function
    s_trans[0,i+1] = y_trans[0,i+1] - c_trans[0,i+1] 

# Plot transition paths
plt.figure(figsize=(9,6))
plt.title('Consumption transition path', fontsize = 20) # title with fontsize 20
plt.plot(c_trans.squeeze(), linewidth=2.0)
plt.ylabel('Consumption')
plt.xlabel('Period')
plt.savefig('c_transpath.png')
plt.show() 

plt.figure(figsize=(9,6))
plt.title('Labor transition path', fontsize = 20) # title with fontsize 20
plt.plot(h_trans.squeeze(), linewidth=2.0)
plt.ylabel('Labor')
plt.xlabel('Period')
plt.savefig('l_transpath.png')
plt.show() 

plt.figure(figsize=(9,6))
plt.title('Output transition path', fontsize = 20) # title with fontsize 20
plt.plot(y_trans.squeeze(), linewidth=2.0)
plt.ylabel('Output')
plt.xlabel('Period')
plt.savefig('o_transpath.png')
plt.show() 

plt.figure(figsize=(9,6))
plt.title('Savings transition path', fontsize = 20) # title with fontsize 20
plt.plot(s_trans.squeeze(), linewidth=2.0)
plt.ylabel('Savings')
plt.xlabel('Period')
plt.savefig('s_transpath.png')
plt.show() 