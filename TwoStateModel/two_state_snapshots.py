# --------------------------------------------------------------------
# INTRODUCTION TO MONTE CARLO WITH PYTHON
# March 2025
# Example program by Jordi Faraudo
# Simulation of a sequence of configurations for an ideal 2 state system
# Calculation of Energy and magnetization of a configuration
# ---------------------------------------------------------------------
#
# Here we import the numpy mathematical library and the plots library
# as in the other examples in the course
#
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

#
# Function implementing the method (Metropolis Monte Carlo) and model
#      The state with s=-1 has zero energy 
#      The state with s=+1 has energy 1
# 
def mcmove(config, N, beta):
        #Loop with a size equal to spins in the system
        for i in range(N):
            for j in range(N):
                    #Select randomly a lattice node
                    a = np.random.randint(0, N)
                    b = np.random.randint(0, N)
                    ## Read state
                    s =  config[a, b]
                    #calculate energy cost of fliping current state
                    if (s<0):
                       #cost of exciting the system  
                       cost=1.0
                    else:
                       #cost of decay
                       cost=-1.0
                    #flip spin or not depending on the cost and its Boltzmann factor
                    ## (acceptance probability is given by Boltzmann factor with beta = 1/kBT)
                    if cost < 0:
                        s = s*(-1)
                    elif rand() < np.exp(-cost*beta):
                        s = s*(-1)
                    config[a, b] = s
        # return the new configuration
        return config

#This function makes an image of the spin configurations
def configPlot(f, config, i, N):
        ''' This modules plts the configuration '''
        X, Y = np.meshgrid(range(N), range(N))
        plt.pcolormesh(X, Y, config, vmin=-1.0, vmax=1.0, cmap='RdBu_r');
        plt.title('MC iteration=%d'%i);
        plt.axis('tight')
        plt.pause(0.1)

#This function calculates the energy of a given configuration for the plots of Energy as a function of T
def calcEnergy(config):
    '''Energy of a given configuration'''
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            State = config[i,j]
            if (State>0):
               energy = energy + 1.0
    return energy

#This function calculates the magnetization of a given configuration
def calcMag(config):
    '''Magnetization of a given configuration'''
    mag = np.sum(config)
    return mag
#
# MAIN PROGRAM
#
#
#  Here we set initial conditions and control the flow of the simulation
#
#size of the lattice
N = 64
#Enter data for the simulation
print("MC Simulation two State system")
print("------------------------------")
print("Epsilon = Energy exited state")
temp = float(input("\n Please enter kT/Epsilon (suggestion 0.5): "))
msrmnt = int(input("\n Enter number of Monte Carlo (Metropolis) iterations (suggestion 100):"))

#Init Magnetization and Energy
step=[]
M=[]
E=[]

#Generate initial condition random state for all sites
config = 2*np.random.randint(2, size=(N,N))-1

#Generate initial condition system in the ground state
#config = np.zeros([N,N]) -1.0

#Calculate initial value of magnetization and Energy
Ene = calcEnergy(config)/(N*N)     # calculate average energy
Mag = calcMag(config)/(N*N)        # calculate average magnetisation
t=0
print('MC step=',t,' Energy=',Ene,' M=',Mag)
#Update 
step.append(t)
E.append(Ene)
M.append(Mag)

#Show initial condition
print('Initial configuration:')
print(config)
#f = plt.figure(figsize=(15, 15), dpi=80);
f = plt.figure(dpi=100)
configPlot(f, config, 0, N)
plt.show()

#Turn on interactive mode for plots
print("Starting MC simulation")
plt.ion()

#Perform the MC iterations
for i in range(msrmnt):
            #call MC calculation
            mcmove(config, N, 1.0/temp)
            #update variables
            t=t+1                              # update MC step
            Ene = calcEnergy(config)/(N*N)     # calculate average energy
            Mag = calcMag(config)/(N*N)        # calculate average magnetisation
            #Update 
            step.append(t)
            E.append(Ene)
            M.append(Mag)

            #plot certain configurations
            if t%10 == 0:
                print('\nMC step=',t,' Energy=',Ene,' M=',Mag)
                print(config)
                configPlot(f, config, t, N)

#Print end
print('\nSimulation finished after',t, 'MC steps')

#interactive plotting off
plt.ioff()

#Show final configuration
configPlot(f, config, t, N)
plt.show()

#Plot evolution of Energy and average state during the simulation
plt.subplot(2, 1, 1)
plt.plot(step, E, 'r+-')
plt.ylabel('Energy')

plt.subplot(2, 1, 2)
plt.plot(step, M, 'b+-')
plt.ylabel('Average State')
plt.xlabel('MC step')

#Show the plot in screen
plt.show()
