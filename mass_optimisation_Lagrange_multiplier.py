import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

def constraint_equation(x, v, I_sp, epsilon):
    # v given in km/s for better convergence of numerical solver
    g0 = 9.81
    c = g0*I_sp/1000
    #because np.sum(c * np.log( (x*c - 1) / (epsilon*x*c) )) is an increasing function 
    #we can set undefined values to zero for numerical approximations using bisection like methods
    #way around using newtons method with crucial initial guesses
    if x <= 0 or np.any((x*c - 1) <= 0):
        return -1
    else:
        return np.sum(c * np.log( (x*c - 1) / (epsilon*x*c) )) - v

def Lagrange_mass_optimisation(v_end, I_sp, epsilon):
    #v_end is given in km/s
    #Ensure arrays are numpy arrays
    I_sp = np.array(I_sp)
    epsilon = np.array(epsilon)
    g0 = 9.81
    c = g0*I_sp/1000 #has units km/s
    #set payload mass to unity
    m_payload = 1
    #safe mass of payload
    m0 = m_payload
    
    
    #Solve the constraint equation
    result = root_scalar(constraint_equation, args = (v_end, I_sp, epsilon), method = 'brentq', bracket = [0.1,1e3])
    if result.converged:
        mu = result.root
    else:
        raise ValueError("Root finding does not converge")
    #print("mu:",mu)
    x0 = mu
    #print("constraint terms",c * np.log( (x0*c - 1) / (epsilon*x0*c) ))
    n = (mu*c - 1) / (epsilon*mu*c)

    #get number of stage
    s = len(I_sp)
    #masses for the seperate stages
    m = np.zeros(s)
    #loop backwards from s - 1 to 0
    for k in range(s-1,-1,-1):
        m[k] = (1- n[k]) / (epsilon[k]*n[k] -1) * m_payload

        #Update payload for lower stage
        m_payload = m_payload + m[k]

    m_empty = epsilon*m
    m_fuel = m - m_empty

    m0 = m0 + np.sum(m)

    return m0, m, m_empty, m_fuel

v_end = 10
I_sp = [300, 300, 300]
epsilon = [0.1, 0.1, 0.1]

#for given payload mass the total mass is given by multiplying by m_payload
m0, m, m_empty, m_fuel = Lagrange_mass_optimisation(v_end, I_sp, epsilon)
print(1/m0*m)

#setting fontsize
fontsize_label = 14
fontsize_legend = 14
fontsize_ticks = 14

epsilon_array = np.linspace(0.05, 0.15, 11)
I_sp_array = [300]

#number of stages
s = 3

for x in I_sp_array:
    I_sp = x*np.ones(s)
    m_normalised = np.zeros((s, len(epsilon_array)))

    for j, eps in enumerate(epsilon_array):
        epsilon = eps*np.ones(s)
        m0, m, m_empty, m_fuel = Lagrange_mass_optimisation(v_end, I_sp, epsilon)

        #normalise m
        m = 1/m0*m 

        #Store values
        m_normalised[:,j] = m

    #Plot the normalised values
    for k in range(s):
        plt.plot(epsilon_array, m_normalised[k, :])
    print(m_normalised)

plt.show()