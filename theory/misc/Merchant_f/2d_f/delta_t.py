import os
import numpy as np
import sys

if __name__ == '__main__':
	time_p = []
	for i in [-6, -5, -4, -3, -2, -1]:
		for j in [1, 2.5, 5, 7.5]:
			time_p.append(j*10**(i))
	print(time_p)
	tau = 0.777660157519
	for idx, val in enumerate(time_p):
		time_p[idx] = val / tau
	print(time_p)
