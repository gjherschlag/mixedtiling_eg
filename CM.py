from   universe import *

#set system parameters
noruns = 20
for i in range(0,noruns):
	maxpp = 50;
	minpp = 0.5;
	dlog  = (np.log(maxpp)-np.log(minpp))/noruns
	pco   = np.exp(dlog*i)*0.5
	system = universe("./parameters.txt",pco=pco, pco_indx=i, tiling=2)
	system.loadProgram("solverscpp")
	system.execute_mixed()