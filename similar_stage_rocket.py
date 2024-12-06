import numpy as np
import matplotlib.pyplot as plt

'''
Similar stage rockets have identical ratios lambda = m_Pl/(m_S + M_F), epsilon = m_S/(m_S + M_F)
and n = m_0/(m_S + m_F) for each stage which follows directly from the first two
Moreover they all have the same specific impulse I_sp
'''

g0 = 9.81    # in m/s
I_sp = 300   # in s

r = 0.05     # ratio of payload/total mass
epsilon = 0.15

num_stages = 10
stages = range(1, num_stages+1)
v_stages = np.array([])

for n in stages:
    v_end = I_sp*g0*n*np.log(1/(r**(1/n) * (1 - epsilon) + epsilon))
    v_stages = np.append(v_stages, v_end)

v_max = I_sp*g0*(1-epsilon)*np.log(1/r)
plt.scatter(stages, v_stages/v_max, label = f'I_sp = {I_sp}, $\epsilon = $ {epsilon}, r = {r}')
plt.xlabel('Number of stages')
plt.ylabel('Final velocity to maximal velocity')
plt.title('Final velocity for stage similar rocket per stage')
plt.legend()
plt.show()