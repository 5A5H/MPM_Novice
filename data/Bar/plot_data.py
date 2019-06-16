#!/anaconda3/envs/PYTHON3/bin/python

# This file is used to plot results
import matplotlib.pyplot as plt

# Load in results
T=[];
VX=[];
VY=[];
VZ=[];
with open('/Users/sash/mpm_2d/data/Bar/MyCVS_V') as file:
    for line in file:
        T.append(   float(line.split(',')[0]) )
        VX.append(  float(line.split(',')[1]) )
        VY.append(  float(line.split(',')[2]) )
        VZ.append(  float(line.split(',')[3]) )

X=[];
Y=[];
Z=[];
with open('/Users/sash/mpm_2d/data/Bar/MyCVS_X') as file:
    for line in file:
        X.append(  float(line.split(',')[1]) )
        Y.append(  float(line.split(',')[2]) )
        Z.append(  float(line.split(',')[3]) )

plt.plot(T, VX, color='blue',label='VX')
plt.plot(T, VY, color='black',label='VY')
plt.plot(T, VZ, color='green',label='VZ')
plt.ylabel('velocity')
plt.xlabel('t-time')
plt.legend(loc='lower right')
plt.show()

plt.plot(T, X, color='blue',label='X')
plt.plot(T, Y, color='black',label='Y')
plt.plot(T, Z, color='green',label='Z')
plt.ylabel('coor')
plt.xlabel('t-time')
plt.legend(loc='lower right')
plt.show()
