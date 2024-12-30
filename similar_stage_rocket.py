import numpy as np
import matplotlib.pyplot as plt

'''
Similar stage rockets have identical ratios lambda = m_Pl/(m_S + M_F), epsilon = m_S/(m_S + M_F)
and n = m_0/(m_S + m_F) for each stage which follows directly from the first two
Moreover they all have the same specific impulse I_sp
'''

g0 = 9.81    # in m/s
I_sp = 300   # in s

r_array = [0.01, 0.05]     # ratio of payload/total mass
epsilon_array = [0.05, 0.1, 0.15]

num_stages = 10
stages = range(1, num_stages+1)



for r in r_array:
    v_differences = np.zeros((len(epsilon_array), num_stages))
    i = 0
    for epsilon in epsilon_array:
        v_stages = np.array([])

        for n in stages:
            v_end = I_sp*g0*n*np.log(1/(r**(1/n) * (1 - epsilon) + epsilon))
            v_stages = np.append(v_stages, v_end)

        v_max = I_sp*g0*(1-epsilon)*np.log(1/r)
        
        #Calculate velocity gain per stage
        v_differences[i] = np.insert(np.diff(v_stages), 0, v_stages[0])
        #Normalise to maximal velocity
        v_differences[i] = v_differences[i]/v_max

        #Update i for next row
        i = i+1


        plt.scatter(stages, v_stages/v_max, label = f'$\epsilon = $ {epsilon}, $\lambda$ = {r}')
        plt.plot(stages, v_stages/v_max)

    plt.xlabel('Number of stages')
    plt.ylabel('Final velocity to maximal velocity')
    plt.title('Final velocity for stage similar rocket')
    plt.legend()
    plt.figure()
    
    
    for k in range(len(epsilon_array)):
        plt.scatter(stages, v_differences[k], label = f'$\epsilon = $ {epsilon_array[k]}, $\lambda$ = {r}')
        #plt.plot(stages, v_differences[k])
    
    plt.xlabel('Number of stages')
    plt.ylabel('Velocity gain to the last stage')
    plt.legend()
    plt.figure()

'''
v_stages = np.array([])

for n in stages:
    v_end = I_sp*g0*n*np.log(1/(r**(1/n) * (1 - epsilon) + epsilon))
    v_stages = np.append(v_stages, v_end)

v_max = I_sp*g0*(1-epsilon)*np.log(1/r)
#Notice the quotient of v_end/v_max is independent of g0 and I_sp but equal I_sp is still assumed
plt.scatter(stages, v_stages/v_max, label = f'$\epsilon = $ {epsilon}, r = {r}')
plt.plot(stages, v_stages/v_max)
plt.xlabel('Number of stages')
plt.ylabel('Final velocity to maximal velocity')
plt.title('Final velocity for stage similar rocket')
plt.legend()
'''

#Calculating necessary payload/total mass ratio for target velocity
v_target = 10*10**3 #target velocity (m/s)
r_array = []

for n in stages:
    r = ((np.exp(-v_target / (I_sp*g0*n) ) - epsilon) / (1 - epsilon))**n
    r_array.append(r)

print(f"{'Stages':<10}{'r_array':<10}")
print("-" * 20)

# Rows
for stage, r in zip(stages, r_array):
    print(f"{stage:<10}{r:<10.3f}")

#Convert arrays to numpy arrays
stages = np.array(stages)
r_array = np.array(r_array)

#Extract negative values of indices (those are not realisable)
r_array = np.where(r_array < 0, np.nan, r_array)
 

#plt.figure()
plt.scatter(stages, r_array, label = fr'$I_{{sp}}$ = {I_sp}, $\epsilon =$ {epsilon}, $v_{{target}}$ = {v_target/1000} km/s')
plt.plot(stages, r_array)
plt.xlim(stages[0], stages[-1])
plt.xlabel('Number of stages')
plt.ylabel('Payload/total mass ratio')
plt.title('r ratio needed for n stages')
plt.legend()

plt.show()