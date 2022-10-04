# 1D linear convection equation: \partial_t u + c * \partial_x u = 0
  
solves linear convection eqn using finite difference methods and periodic bc.
Initial condition is given by Gaussian function.
amplitude * \exp( - alpha * (x-x_0)^2 )

PETSc DMDA is used for vector parallelization and TS for time stepping.

To run:
```
make
./exe
```
run by specifying number of points
```
make
./exe -da_grid_x 500
```



To run in parallel
```
mpirun -np 2 ./exe
```
Then plot solution in gnuplot
```
gnuplot>plot 'solution_0000000.dat' u 1:2 w l lw 4, 'solution_0000050.dat' u 1:2 w l lw 4
```
