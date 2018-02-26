# FSIM
These are codes for simulations and real data analysis, which appears in the manscript entitled, "Functional Single Index Model" by Jiang, F., Baek, S., Cao, J., and Ma, Y. 
All files ended with .f90 is Fortran code, which should be compiled first. 

In the folder "organized", 
cvd1.dat -- reponse variable 
x1.data -- covariates which are correspondent to one-year lagged response
simu1-beta.f90 -- estimation for Simulation 1 for beta
simu1-gamma.f90 -- estimation for Simulation 1 for gamma
simu2-beta.f90 -- estimation for Simulation 2
simu3-beta.f90 -- estimation for Simulation 3 for beta
simu3-gamma.f90 -- estimation for Simulation 3 for gamma
real.f90 -- estimation for real data

In the folder "stacking", there are computational codes for stacking method, which is described in section 4 in the manuscript. 
simu1c.f90 -- estimation for Simulation 1
simu2c.f90 -- estimation for Simulation 2
simu3c.f90 -- estimation for Simulation 3
simu1c.txt -- estimation results of Simulation 1
simu2c.txt -- estimation results of Simulation 1
simu3c.txt -- estimation results of Simulation 1


