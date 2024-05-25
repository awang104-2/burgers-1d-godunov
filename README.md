# Burgers 1D FV Solver (Godunov Method)

Solves the Burgers 1D equation using the Godunov method with parameters written in the code. To use the Godunov method, the program discretizes the continuous x-range it's solving over into small, discrete intervals, then solves a localized Riemann problem for the inviscid Burgers equation at each interval, with initial conditions written in the python file. The first few graphs in the tester demonstrate the Riemann problem solver, and the last graph shows the solutions to the Burgers 1D equation using the Godunov method.
 
