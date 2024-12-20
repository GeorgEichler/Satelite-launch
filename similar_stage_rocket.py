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
#Notice the quotient of v_end/v_max is independent of g0 and I_sp but equal I_sp is still assumed
plt.scatter(stages, v_stages/v_max, label = f'$\epsilon = $ {epsilon}, r = {r}')
plt.xlabel('Number of stages')
plt.ylabel('Final velocity to maximal velocity')
plt.title('Final velocity for stage similar rocket per stage')
plt.legend()


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

plt.figure()
plt.scatter(stages, r_array, label = fr'$I_{{sp}}$ = {I_sp}, $\epsilon =$ {epsilon}, $v_{{target}}$ = {v_target/1000} km/s')
plt.xlabel('Number of stages')
plt.ylabel('Payload/total mass ratio')
plt.title('r ratio needed for n stages')
plt.legend()

plt.show()