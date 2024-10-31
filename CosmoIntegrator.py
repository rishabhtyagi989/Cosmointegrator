#this code solves for the scale factor equation 
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#defining all the relevant constants in natural units i.e hbar=c=kb=1

H_0 = 0.6688*(1/3.0856)*10**(-17)
om_m_0 = 0.321
om_r_0 = 0.321/3412
om_l_0 = 1-om_m_0-om_r_0

#Solves the Scale factor equation for real universe      
def SFET(Y, t):
    a = Y[0]
    w = Y[1]
    dadt = -1*w
    dwdt= (H_0**2)*(om_r_0*a**(-4)+(1/2)*om_m_0*a**(-3)-1*om_l_0)*a
    return [dadt, dwdt]


# Initial conditions
Y_0 = [1,H_0]

# time points
t = np.linspace(0, 13.84e9, num=10000)*(3600*24*365)

#Solve the scale factor equations
Ysol = odeint(SFET,Y_0,t)
a = Ysol[:, 0]

#plot scale factor vs time
plt.plot(t/(3600*24*365),a,'g-',linewidth=2,label='Rad+Mat+DE')
plt.ylim(-0.1,1.1)
plt.axhline(0)
plt.xlabel('time (Yrs)')
plt.ylabel('$a(t)$')
plt.legend()
plt.show();
