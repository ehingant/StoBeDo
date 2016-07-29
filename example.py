#!/usr/bin/python
from stobedo import *
import matplotlib.pyplot as plt
from matplotlib import rc


# First exemple: Nucleation time 
# T_N = inf { t > 0 : C_N(t) = 1 }
# N: size of the nucleus	
# Xp: number of paths (simulations)
# M: Number of particles
# iniial condition C_1(0) = M, C_i(0)=0 for all i >1, surely.

M = 1000
Xp = 100
N = 100
T = np.zeros(Xp)

for n in range(0,Xp):
	T[n] = CalculNucleationTime(M,N)

plt.figure(1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.title('Nucleation time $T_N =\inf\{ t\geq 0 : C_N(t)=1\}$ for '+str(Xp)+'simulations')
plt.ylabel('Normalized number')
plt.xlabel('$T_N$')
plt.hist(T,50,normed=1,facecolor='green',alpha=0.75)

#Second exemple: Sample path of C_i for some i on [0,Tmax)
M = 1000
Xp = 10
i = 3
Tmax = 1

plt.figure(2)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.title('Sample path of $C_i$ with i='+str(i)+' for '+str(Xp)+'  simulations')
plt.ylabel('Normalized $C_i(t)$')
plt.xlabel('time $t$')
for n in range(0,Xp):
	[x,y] = CalculPath(M,i,Tmax)
	plt.plot(x,y,drawstyle='steps')

#Third exemple: Distribution of the C_i's  at time T

M = 1000
Xp = 10
i = 3
T = 1

plt.figure(3)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.title('Distirbution of $C$  at time t='+str(T)+' for '+str(Xp)+'  simulations')
plt.ylabel('Number of clusters')
plt.xlabel('size $i$')
for n in range(0,Xp):
	X,Y = CalculDistribution(M,T)
	plt.plot(X,Y)

plt.show()

