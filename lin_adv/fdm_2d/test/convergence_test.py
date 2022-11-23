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
lc = [40, 60, 80]
for n in range(len(lc)):
   print('Running with Np = '+str(lc[n]))
   os.system(exe+' -da_grid_x '+str(lc[n])+' -da_grid_y '+str(lc[n])+' -ts_type rk -ts_rk_type 4 > out')
   os.system('tail -1 out >> error.dat')

# read error from file
data = np.loadtxt('error.dat')
plt.autoscale(tight=False)
dx = data[:,0]
err_l2 = data[:,1]
plt.loglog(dx, err_l2,'o-', linewidth=2)
plt.loglog(dx, dx, '-', linewidth=2)
plt.loglog(dx, dx*dx,'-', linewidth=2)
plt.loglog(dx, dx*dx*dx,'-', linewidth=2)
plt.loglog(dx, dx*dx*dx*dx,'-', linewidth=2)
plt.xlabel('dx')
plt.ylabel('Error')
plt.legend(('L2', '1st', '2nd', '3rd', '4th'))
plt.savefig('error.pdf')
plt.show()
