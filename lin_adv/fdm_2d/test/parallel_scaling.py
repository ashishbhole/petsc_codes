import os
import numpy as np
from matplotlib import pyplot as plt
import argparse

os.system('sh ./speedup.sh "1 2 3 4 5 6 7 8 9 10"') 

# read error from file
data = np.loadtxt('times.txt')
plt.autoscale(tight=False)
np = data[:,0]
time = data[:,1]
plt.plot(np, time,'o-', linewidth=2)
plt.xlabel('No of processors')
plt.ylabel('Time')
plt.savefig('time.pdf')
plt.show()
