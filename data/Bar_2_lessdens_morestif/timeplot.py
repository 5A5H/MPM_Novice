
# This python script simply generates a plot out of a two column cvs file

import sys
import matplotlib.pyplot as plt

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

# Load in results
T=[];
VX=[];
with open(sys.argv[1]) as file:
    for line in file:
        T.append(   float(line.split(',')[0]) )
        VX.append(  float(line.split(',')[1]) )

plt.plot(T, VX)
plt.ylabel('Y')
plt.xlabel('t-time')
plt.show()
