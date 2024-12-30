import numpy as np
import matplotlib.pyplot as plt

'''
Similar stage rockets have identical ratios lambda = m_Pl/(m_S + M_F), epsilon = m_S/(m_S + M_F)
and n = m_0/(m_S + m_F) for each stage which follows directly from the first two
Moreover they all have the same specific impulse I_sp
'''


# Set font sizes for plots
font_size_labels = 15
font_size_ticks = 15
font_size_legend = 15

#Change figuresize for plots globally
plt.rcParams['figure.figsize'] = [8, 6]

g0 = 9.81    # in m/s
I_sp = 300   # in s

r_array = [0.01, 0.05]     # ratio of payload/total mass
epsilon_array = [0.05, 0.1, 0.15]

num_stages = 8
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

    plt.xlabel('Number of stages', fontsize = font_size_labels)
    plt.xticks(fontsize = font_size_ticks)
    plt.ylabel('Final velocity/maximal velocity', fontsize = font_size_labels)
    plt.yticks(fontsize = font_size_ticks)
    #plt.title('Final velocity for stage similar rocket')
    plt.legend(fontsize = font_size_legend)
    plt.figure()
    
    
    for k in range(len(epsilon_array)):
        plt.scatter(stages[1:], v_differences[k][1:], label = f'$\epsilon = $ {epsilon_array[k]}, $\lambda$ = {r}')
        #plt.plot(stages, v_differences[k])
    
    plt.xlabel('Number of stages', fontsize = font_size_labels)
    plt.xticks(fontsize = font_size_ticks)
    plt.ylabel('Velocity ratio gain compared to previous stage', fontsize = font_size_labels)
    plt.yticks(fontsize = font_size_ticks)
    #plt.ylim([0,0.3])
    plt.legend(fontsize = font_size_legend)
    plt.figure()




#Calculating necessary payload/total mass ratio for target velocity
v_target = 9*10**3 #target velocity (m/s)
I_sp_array = [300]   # in s

stages = np.array(stages)

for I_sp in I_sp_array:
    lambda_differences = np.zeros((len(epsilon_array), num_stages))
    i = 0
    for epsilon in epsilon_array:
        lambda_array = np.array([])

        for n in stages:
            lambda_end = ((np.exp(-v_target / (I_sp*g0*n) ) - epsilon) / (1 - epsilon))**n
            lambda_array = np.append(lambda_array, lambda_end)

        lambda_array = np.where(lambda_array < 0, np.nan, lambda_array)

        #Calculate payload ratio gain per stage
        lambda_differences[i] = np.insert(np.diff(lambda_array), 0, lambda_array[0])

        i = i+1

        plt.scatter(stages, lambda_array, label = fr'$I_{{sp}}$ = {I_sp}, $\epsilon =$ {epsilon}')
        plt.plot(stages, lambda_array)

    plt.xlabel('Number of stages', fontsize = font_size_labels)
    plt.xticks(fontsize = font_size_ticks)
    plt.ylabel('Payload to total mass ratio $\lambda$', fontsize = font_size_labels)
    plt.yticks(fontsize = font_size_ticks)
    #plt.title('$\lambda$ ratio needed for n stages')
    plt.legend(fontsize = font_size_legend)
    plt.figure()

    for k in range(len(epsilon_array)):
        plt.scatter(stages[1:], lambda_differences[k][1:],
                    label = fr'$I_{{sp}}$ = {I_sp}, $\epsilon =$ {epsilon_array[k]}')
        plt.plot(stages[1:], lambda_differences[k][1:])
    
    plt.xlabel('Number of stages', fontsize = font_size_labels)
    plt.xticks(fontsize = font_size_ticks)
    plt.ylabel('$\lambda$ gain compared to previous stage', fontsize = font_size_labels)
    plt.yticks(fontsize = font_size_ticks)
    #plt.xlim(0, num_stages)
    #plt.ylim([0,0.3])
    plt.legend(fontsize = font_size_legend)
    plt.figure()

#closes all unused figures
plt.close()
plt.show()
exit()

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