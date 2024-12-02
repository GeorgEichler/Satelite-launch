import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def vertical_rocket_model(t, y, m0, ve, dmdt, Cd, R):
    #Input variables:
    #t    - time
    #y    - state vector (z,v)
    #m0   - initial mass of rocket
    #ve   - velocity of expelled gas
    #dmdt - speed of mass loss
    #Cd   - drag coefficient
    #R    - rocket radius
    
    # Constants
    G = 6.674*10**(-11) #gravitational constant
    Re = 6.378*10**6    #radius of the earth
    M = 5.972*10**24    # mass of the earth
    rho0 = 1.2*10**3    #air density on sea level
    H = 8.333*10**3     #height scale
    A = np.pi*R**2      #cross-sectional area of the rocket

    # Unpack the state variables
    z, v = y
    dzdt = v 
    dvdt = ve/(m0 - dmdt*t)*dmdt -G*M/(Re + z)**2 - 1/(2*(m0 - dmdt*t))*Cd*A*rho0*np.exp(-z/H)*v**2 

    return [dzdt, dvdt]

# Initial conditions

m0 = 137*10**3 + 2169*10**3 + 141*10**3   # empty mass + fuel mass + payload (data taken from Wikipedia Saturn V)
m_final = 137*10**3 + 141*10**3           #final mass
ve = 3*10**3                              # velocity of expelled gas
dmdt = 12.5*10**3                         # speed of mass loss (kg/s)
Cd = 0                                    # drag coefficient
R = 10                                    # rocket radius (m)

#Calculate burn time
tb = (m0 - m_final)/dmdt

print(f'The burning time is {tb}s')

# Initial state vector
y0 = [0, 0]

# Time interval
t_span = [0, tb]
t_eval = np.linspace(t_span[0], t_span[1], 1000)

sol = solve_ivp(vertical_rocket_model, t_span, y0, t_eval = t_eval, args = (m0, ve, dmdt, Cd, R))

plt.plot(sol.t, sol.y[0], label = 'Rocket height')
plt.xlabel('Time t (s)')
plt.ylabel('Height z (m)')
plt.legend()

plt.figure()

plt.plot(sol.t, sol.y[1], label = 'Rocket velocity')
plt.xlabel('Time t (s)')
plt.ylabel('Velocity v (m/s)')
plt.legend()
plt.show()
