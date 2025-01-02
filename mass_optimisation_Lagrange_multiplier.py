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

v_end = 9
I_sp = [300, 300, 300]
epsilon = [0.1, 0.1, 0.1]

#for given payload mass the total mass is given by multiplying by m_payload
m0, m, m_empty, m_fuel = Lagrange_mass_optimisation(v_end, I_sp, epsilon)
print(m0)

#setting fontsize
fontsize_label = 14
fontsize_legend = 14
fontsize_ticks = 14

epsilon_array = np.linspace(0.05, 0.15, 101)
I_sp_array1 = [[250,250],[300,300],[350,350]]
I_sp_array2 = [[250,250,250],[300,300,300],[350,350,350]]

I_sp_array1 = [[300,350]]
I_sp_array2 = [[300,350,350]]

#number of stages
s = 2

m02 = []

for Isp in I_sp_array1:
    I_sp = np.array(Isp)
    m_normalised = np.zeros((s, len(epsilon_array)))

    for j, eps in enumerate(epsilon_array):
        epsilon = eps*np.ones(s)
        m0, m, m_empty, m_fuel = Lagrange_mass_optimisation(v_end, I_sp, epsilon)

        m02.append(m0)

        #normalise m
        m = 1/m0*m 

        #Store values
        m_normalised[:,j] = m

    #Plot the normalised values
    for k in range(s):
        plt.plot(epsilon_array, m_normalised[k, :], label = fr'Stage {k+1}')

    plt.xlabel(fr'Structural ratio $\epsilon$', fontsize = fontsize_label)
    plt.xticks(fontsize = fontsize_ticks)
    plt.ylabel(fr'Stage masses normalised to $m_0$', fontsize = fontsize_label)
    plt.yticks(fontsize = fontsize_ticks)
    plt.legend(fontsize = fontsize_legend)
    #plt.title(f'Isp = {Isp}')
    plt.figure()

#number of stages
s = 3
m03 = []

for Isp in I_sp_array2:
    I_sp = np.array(Isp)
    m_normalised = np.zeros((s, len(epsilon_array)))

    for j, eps in enumerate(epsilon_array):
        epsilon = eps*np.ones(s)
        m0, m, m_empty, m_fuel = Lagrange_mass_optimisation(v_end, I_sp, epsilon)

        m03.append(m0)

        #normalise m
        m = 1/m0*m 

        #Store values
        m_normalised[:,j] = m

    #Plot the normalised values
    for k in range(s):
        plt.plot(epsilon_array, m_normalised[k, :], label = fr'Stage {k+1}')

    plt.xlabel(fr'Structural ratio $\epsilon$', fontsize = fontsize_label)
    plt.xticks(fontsize = fontsize_ticks)
    plt.ylabel(fr'Stage masses normalised to $m_0$', fontsize = fontsize_label)
    plt.yticks(fontsize = fontsize_ticks)
    plt.legend(fontsize = fontsize_legend)
    #plt.title(f'Isp = {Isp}')
    plt.figure()

#Transform to numpy arrays to allow vector operations
m02 = np.array(m02)
m03 = np.array(m03)

mass_ratio = m03[:101]/m02[:101]

plt.plot(epsilon_array, mass_ratio)
plt.xlabel(fr'Structural ratio $\epsilon$', fontsize = fontsize_label)
plt.xticks(fontsize = fontsize_ticks)
plt.ylabel(fr'Mass ratio $m_0^{{(3)}}/m_0^{{(2)}}$', fontsize = fontsize_label)
plt.yticks(fontsize = fontsize_ticks)
#print(m02[:101])
#print(m03[:101])

print(fr'Stage 2: mass for $\epsilon = 0.05$ is {m02[0]} and $\epsilon = 0.15$ is {m02[100]}')
print(f'Ratio for stage 2 of these masses is {m02[100]/m02[0]}')
print(fr'Stage 3: mass for $\epsilon = 0.05$ is {m03[0]} and $\epsilon = 0.15$ is {m03[100]}')
print(f'Ratio for stage 3 of these masses is {m03[100]/m03[0]}')

plt.show()