import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('points.txt', dtype=float, delimiter=', ')
print(data)

# points = data[125:126, :]
#print(data[:,2])
len = data[:,2]	
last = 0
index = 0

# for i in len:
# 	if(i > 125):
# 		if(i == last):
			
# 		else:
# 			last = i #new time step
# 	print(i)
# 	i++