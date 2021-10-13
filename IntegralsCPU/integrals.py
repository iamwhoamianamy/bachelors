from scipy.integrate import quad, dblquad, tplquad
import numpy as np
import math

f3 = lambda x, y, z: x*x + y*y
f2 = lambda x, y: x*x + y*y


#print(tplquad(f3, -1, 1, -1, 1,-1, 1)[0])
print(dblquad(f2, -1, 1, -1, 1)[0])