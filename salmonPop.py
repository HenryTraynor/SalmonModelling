from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# outside constants
C = 6e7 # carrying capacity

def odes(x, t):
    # constants
    a1 = 800          # salmon birth rate
    a2 = .002*C        # rivalry rate
    a3 = 1/365      # adult maturation rate
    a4 = .005          # adult death rate
    a5 = 2.71828**(-((t%365)-238)**2/200)         # sexual maturation rate
    a6 = .001         # harvesting rate
    a7 = 1/66      # spawn death rate
    a8 = 1.00167e-9*C  # salmon-bear predation
    a9 = 1500/C        # number of bears

    # assign each ODE to a vector element
    J = x[0]
    A = x[1]
    S = x[2]

    # define each ODE as follows
    dJdt = a1*S-a2*J**2-a3*J
    dAdt = a3*J-a4*A-a5*A-a6*A
    dSdT = a5*A-a7*S-a8*S*a9

    return [dJdt, dAdt, dSdT]

# initial conditions
x0 = [.75, .125, .125]

# test the defined ODEs
# print(odes(x=x0,t=0))

# declare a time vector (time window)
t = np.linspace(0,3000,1000)
x = odeint(odes,x0,t)

J = x[:,0]
A = x[:,1]
S = x[:,2]

# plot results
plt.semilogy(t,J/(J+A+S))
plt.semilogy(t,A/(J+A+S))
plt.semilogy(t,S/(J+A+S))
plt.gca().legend(('Juveniles','Adults','Spawning'))
plt.xlabel("Time (days)")
plt.ylabel("Population Proportion")
plt.show()
