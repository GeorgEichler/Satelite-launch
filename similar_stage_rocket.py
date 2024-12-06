'''
Similar stage rockets have identical ratios lambda = m_Pl/(m_S + M_F), epsilon = m_S/(m_S + M_F)
and n = m_0/(m_S + m_F) for each stage which follows directly from the first two
Moreover they all have the same specific impulse I_sp
'''

g0 = 9.81    # in m/s
I_sp = 300   # in s
m_p = 10000 #payload mass
m0 = 100000  #total mass

r = m_p/m0
epsilon = 0.1