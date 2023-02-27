# Ising-model
This is an example of Monte Carlo simulation of the Ising model in python.

We have programs for generating configurations by Monte Carlo and for the full simulation of the phase transition, calculating magnitudes as a function of temperature.

In the "snapshots" program we generate a number of Monte Carlo configurations and plot them at the screen (close each window to generate the next configuration)
The ising program performs the full simulation and may require some time to run depending on machine and the number of temperatures sampled in the calculation.

The program employs the numpy and matplotlib libraries.
The code employs python3. The older python2  version is made available as it is in a dedicated folder.

The algorithm is based on Dr Rajesh Singh blog with notebook resources in Physics. 
https://rajeshrinet.github.io/blog/2014/ising-model/
Adapted by Jordi Faraudo for teaching purposes

