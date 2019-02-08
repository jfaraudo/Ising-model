# -----------------------------------------------------------------
# INTRODUCTION TO PYTHON2 SESSION
# May 2018 
# Example 4 program by Jordi Faraudo
# Simulation of a sequence of configurations for the 2D Ising model
# -----------------------------------------------------------------
# 
# The algorithm is based on Dr Rajesh Singh (Cambridge University) blog with python resources in Physics. 
# https://rajeshrinet.github.io/blog/2014/ising-model/
# Adapted by Jordi Faraudo 2018 for teaching purposes
#
# Here we import the numpy mathematical library and the plots library
# as in the other examples in the course
#
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

      
#
# Here we define the interactions of the model (2D spin Ising model) 
# and the solution method (Metropolis Monte Carlo) 
def mcmove(config, N, beta):
        #Loop with a size equal to spins in the system 
        for i in range(N):
            for j in range(N):
                    #generate integer random number between 0 and N            
                    a = np.random.randint(0, N)
                    b = np.random.randint(0, N) 
                    ## Perform a change in the system according to monte carlo move rule
                    s =  config[a, b]
                    #calculate energy cost of this new configuration (the % is for calculation of periodic boundary condition)
                    nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
                    cost = 2*s*nb
                    #flip spin or not depending on the cost and its Boltzmann factor
                    ## (acceptance probability is given by Boltzmann factor with beta = 1/kBT)
                    if cost < 0:	
                        s = s*(-1)
                    elif rand() < np.exp(-cost*beta):
                        s = s*(-1)
                    config[a, b] = s
        # return the new configuration
        return config

    #This function makes the plot                    
                    
def configPlot(f, config, i, N, n_):
        ''' This modules plts the configuration '''
        X, Y = np.meshgrid(range(N), range(N))
        plt.pcolormesh(X, Y, config, vmin=-1.0, vmax=1.0, cmap='RdBu_r');
        plt.title('MC iteration=%d'%i); 
        plt.axis('tight') 
	plt.show()
#
# MAIN PROGRAM
#
#
#  Here we set initial conditions and control the flow of the simulation
#         
#size of the lattice
N = 64 
#Enter data for the simulation
temp = float(raw_input("\n Please enter temperature in reduced units (suggestion 0.4): "))
msrmnt = int(raw_input("\n Enter number of Monte Carlo iterations (suggestion 1000):"))
#Generate initial condition
config = 2*np.random.randint(2, size=(N,N))-1
#Show and Plot initial condition
print 'Initial configuration:'
print config
#f = plt.figure(figsize=(15, 15), dpi=80); 
f = plt.figure(dpi=100);    
configPlot(f, config, 0, N, 1);
     
#Perform the MC iterations
for i in range(msrmnt+1):
            #call MC calculation
            mcmove(config, N, 1.0/temp)
            #plot certain configurations
            print 'MC step=',i
            if i == 1:        
                print config
		configPlot(f, config, i, N, 2)
            if i == 10:
                print config
                configPlot(f, config, i, N, 3)
            if i == 100:      
                print config
                configPlot(f, config, i, N, 4)
            if i == 500:     
                print config
                configPlot(f, config, i, N, 5)
            if i == msrmnt:    
                print config
                configPlot(f, config, i, N, 6)
 
