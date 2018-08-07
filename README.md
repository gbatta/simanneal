# simanneal
Simulated annealing for Julia

This file implements simulated annealing withing Julia. 

The command stucture is as follows:
`simanneal(obj_fn,obj_fn_args::Tuple,lb,ub,c,nn_g_max,nn_a_max,del,max_eval,func_tol,stall_iter_lim_sc=5,rand_init=500990)`

## Arguments

`obj_fn`: The name of the program that returns the cost you are trying to minimize.

`obj_fn_arguments`: Tuple of arguments to the objective function. The first argument must be a vector of initial parameters you are going to test.

`lb`: Vector-valued lower bound of paramater values to search over.

`ub`: Vector-valued upper bound of paramater values to search over.

`c`: Control parameter. Usually advised to range between 1.0 and 10.0.

`nn_g_max`: The number of points to generate before annealing commences.

`nn_a_max`: The number of annealing cycles until reannealing occurs.

`del`: For reannealing, a vector of perturbation sizes for each parameter. A decent way to generate this is to set `del` equal to `abs(ub-lb)*1e-2`

`max_eval`: The maximum number of generated points before the algorithm stops.

`func_tol_array_size` reannealing cycles should fall before we declare the algorithm has converged.

`func_tol_array_size`: The number of reannealing cycles' minimum costs to look back on, to determine convergence.

`noisy`: Set equal to 1 (0 otherwise) if you want each iteration to be outputted to a .csv file entitled "simanneal_output.csv". 

`stall_iter_lim` (OPTIONAL): Minimum number of reannealing cycles.

`rand_init` (OPTIONAL): Integer seed value for random number generators.

`a=(N_cash,c_b,false,actual_moments,raw_weights,moment_number,num_obs,true,true,true,true,true)`

`simanneal(dynamic_programming,a,lb_cash,ub_cash,1,100,10,del_cash,2500,1.0e-3)`
