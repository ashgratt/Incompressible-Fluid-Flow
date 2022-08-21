import numpy as np

import matplotlib.pyplot as plt

# Colebrook-White Equation Parameters

rho = 1000 # [kg / m3]
mu = 0.001 # [Pa s]
d = (50)*10**-3 # [m]
epsilon = (0.025)*10**-3 # [m]
u = 1 # [m / s]

# Reynolds Number Calculation Function:
def Reynolds_Number(rho, u, d, mu):
    
    Re = (rho * u * d)/(mu)
    
    return Re

# Colebrook-White Calculation Function:
def CW(fD, epsilon, Re):
    
    LHS = (1/(fD**0.5))
    
    RHS = 2 * np.log10(   ((epsilon)/(3.7*d))     +     ((2.51)/(Re*(fD**0.5)))   )
    
    Func = LHS + RHS
    
    return Func

# Solver Initialisation:
Re = Reynolds_Number(rho, u, d, mu) 

if Re < 4000:
    
    print('Warning, Reynolds number below validity value of < 4000 (i.e. Re = ' + str(Re) + '. Consider using different correlation.')

fD0 = 0.05
fD1 = fD0*0.99
fD =[fD0, fD1]

loop_counter_array = [0, 1]

i = 1

# Set Tolerance:
tol = 10**-6


while abs(fD[i] - fD[i-1]) > tol:
    
    fDi = fD[i] - CW(fD[i], epsilon, Re) * ((fD[i] - fD[i-1])/(CW(fD[i], epsilon, Re) - CW(fD[i-1], epsilon, Re)))
    
    fD.append(fDi)
    
    i +=1
    
    loop_counter = loop_counter_array[-1] + 1
    
    loop_counter_array.append(loop_counter)


# Show result and provide report on iterations.
print("Result is fD = " + str(fD[-1]))
print("After " + str(loop_counter) + " iterations.")

#print(len(fD))
#print(len(loop_counter_array))

# Plot convergence.
plt.plot(loop_counter_array, fD)
plt.xlabel('Iteration Number')
plt.ylabel('fD')


# Program Critique:

# Secant method is sensitive to poor initial guesses, so consider how one might solve this.
# Obviously no user input information, it's all defined within the py file - consider user input prompt or interface
# Similar to PEL, along with roughness, allow for drop-down selection of roughnesses etc.
# Warnings to tell user for Reynolds number which CW equation is not valid, better yet, automatic detection and switching of correlation.