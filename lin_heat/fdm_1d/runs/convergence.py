import os
import numpy as np
from matplotlib import pyplot as plt
import argparse

# get current directory 
current_dir = os.getcwd()

# path to executable
exe = '../src/exe'
if os.path.isfile(exe)==False:
    print("Could not find ", exe)
    exit()

# no of points
command='rm error.dat'
os.system(command)
lc = [50, 100, 200, 400, 800]
for n in range(len(lc)):
   print('Running with Np = '+str(lc[n]))
   os.system('mpirun -n 2 '+exe+' -da_grid_x '+str(lc[n])+' -ts_type ssp > out')
   os.system('tail -1 out >> error.dat')

# read error from file
data = np.loadtxt('error.dat')
plt.autoscale(tight=False)
dx = 1/data[:,0]
err_l2 = data[:,1]
plt.loglog(dx, err_l2,'o-', linewidth=2)
plt.loglog(dx, dx,'-', linewidth=2)
plt.loglog(dx, dx*dx,'-', linewidth=2)
plt.loglog(dx, dx*dx*dx,'-', linewidth=2)
plt.loglog(dx, dx*dx*dx*dx,'-', linewidth=2)
plt.xlabel('dx')
plt.ylabel('Error')
plt.legend(('L2', '1st', '2nd', '3rd', '4th'))
plt.savefig('error.pdf')
plt.show()
