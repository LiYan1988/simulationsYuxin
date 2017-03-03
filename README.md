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
	