1. When link lengths are very long, TR and GN should not exchange their MIP start.
	This is because there are some links longer than TR but shorter than GN 
	transmission reach. So routings given by GN and TR are different.
	But when the link length are small, it's beneficial to exchange MIP start 
	between TR and GN algorithms.
2. The storage of Gurobi model instances may use a lot of memory, discard them to
	save memory.
3. Given variables a general lb/ub, do not make the Nmax mistakes
4. Always be careful about the numerical issues. It's necessary to post process
	solution values and avoid long trailing decimals. All the numbers should be 
	rounded to a fixed accuracy.
5. Always do small experiments on cluster before large simulations to find bugs and
	save resources.
6. Interior-point method may lead to worse numerical results, due to the long 
	trailing effect of the solutions.
7. Cluster testing procedure:
	a. test python/matlab/C++ file on local machine with a small instance
	b. test bash/slurm file on cluster with one or multiple small instances
	c. write a prepare file that produces all simulation files and a python/bash
		that run all these files at once, if necessary
	d. finally test the run_batch file
8. To control the memory of the allocated nodes, use --mem in sbatch, for example,
	--mem=64G. This will allocate a node with memory larger than the specified 
	value.
9. To compare with Li's algorithm, only one modulation format is needed. And we 
	should set Nmax=0 in the CLGN. This is because Li's MILP does not consider 
	regenerater.
10. How to check actual CPUTime on Rivanna?
	- sacct -S <starttime> -E <endtime> --format=CPUTimeRAW > cputime.txt
	- then download cputime.txt, and use matlab to calculate the amount of time 
	- it seems that the economy discount rate is 0.125
11. slurm sbatch parameters:
	- --cpus-per-task and --ntasks-per-node together control how many cpus are 
		allocated. In default, --cpus-per-task is 1, so --ntasks-per-node will control
		the number of cpus allocated. Otherwise, the number of allocated cpu equals to
		the product of --cpus-per-task and --ntasks-per-node.
	- --mem-per-cpu and --mem control the memory allocated. But they are mutual 
		exclusive
	- Rivanna charges according to max{cpus, memory/6G}. So it's impossible to
		fool the system to charge less than what you use
12. It is possible to directly modify the model to change the bounds of variables,
	  and thus update demands_added and demands_fixed without rewrite the model again.
	  Modifying the model can also automatically provide the solution from the previous
	  solving as the MIP start, which can avoid bugs (if numerical stability is not a 
	  problem).
13. In the gurobi.log file, I found that when GN is severely worse than TR, the GN does 
	  not accept the MIP start given by TR. So change MIP start so that TR gives TR and 
	  GN gives GN.
14. Set presolve to 2?
	  Yes
15. squeue commonly used command:
	  squeue -u ly6j --Format=jobid:.7,name:.14,partition:.11,state:.8,starttime:.20,timeused:.12,endtime:.20,reasonlist:.17,numcpus:.3,minmemory:.11
