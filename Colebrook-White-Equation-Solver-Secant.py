# Colebrook-White Equation Solver - Secant Method

# This file takes in-code inputs of density (rho), viscosity (mu), inner diameter (d), absolute pipe roughness (epsilon),
# average velocity and an initial guess of Friction Factor (fD0).

# It then calculates the Darcy Friction Factor (fD), rounded to 5 decimal places, using the Colebrook-White Equation, which it solves
# via the Secant method. 

# The file reports the Darcy Friction Facotr in the console, states how many iterations were
# required to converge, and then plots the convergence of fD. 

# If the file does not converge, it is suggested that the Friction Factor Estimate is amended - if too large, it will
# not converge easily. A generally good estimate is 0.05 - Friction Factors are generally in the range of 0.001 to 0.1,
# according to most Moody diagrams.

# Program Critique:

# Secant method is sensitive to poor initial guesses, so consider how one might solve this.
# Obviously no user input information, it's all defined within the py file - consider user input prompt or interface
# Similar to PEL, along with roughness, allow for drop-down selection of roughnesses etc.
# Warnings to tell user for Reynolds number which CW equation is not valid, better yet, automatic detection and switching of correlation.

import numpy as np

import matplotlib.pyplot as plt

# Colebrook-White Equation Parameters

rho = 1000 # [kg / m3]
mu = 0.001 # [Pa s]
d = (50)*10**-3 # [m]
epsilon = (0.025)*10**-3 # [m]
u = 1 # [m / s]

# Set Initial Guess at Friction Factor:
fD0 = 0.05 # [-]

# Set Tolerance:
tol = 10**-6

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

fD1 = fD0 * 0.999
fD =[fD0, fD1]

loop_counter_array = [0, 1]

i = 1

while abs(fD[i] - fD[i-1]) > tol:
    
    fDi = fD[i] - CW(fD[i], epsilon, Re) * ((fD[i] - fD[i-1])/(CW(fD[i], epsilon, Re) - CW(fD[i-1], epsilon, Re)))
    
    fD.append(fDi)
    
    i +=1
    
    loop_counter = loop_counter_array[-1] + 1
    
    loop_counter_array.append(loop_counter)


# Show result and provide report on iterations.

Result = round(fD[-1], 5)

print("Result is fD = " + str(Result))
print("After " + str(loop_counter) + " iterations.")

#print(len(fD))
#print(len(loop_counter_array))

# Plot convergence.
plt.plot(loop_counter_array, fD)
plt.xlabel('Iteration Number')
plt.ylabel('fD')
