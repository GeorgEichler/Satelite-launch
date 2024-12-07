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

def Lagrange_mass_optimisation(v_end, I_sp, epsilon, m_payload):
    #Ensure arrays are numpy arrays
    I_sp = np.array(I_sp)
    epsilon = np.array(epsilon)
    g0 = 9.81
    c = g0*I_sp/1000 #has units km/s
    #safe mass of payload
    m0 = m_payload
    #x0 = 0.1
    
    #Solve the constraint equation
    result = root_scalar(constraint_equation, args = (v_end, I_sp, epsilon), method = 'brentq', bracket = [0.1,1e3])
    if result.converged:
        mu = result.root
    else:
        raise ValueError("Root finding does not converge")
    
    n = (1 + mu*c) / (epsilon*mu*c)

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
I_sp = [400, 350, 300]
epsilon = [0.1, 0.15, 0.2]
m_payload = 5000

m0, m, m_empty, m_fuel = Lagrange_mass_optimisation(v_end, I_sp, epsilon, m_payload)
