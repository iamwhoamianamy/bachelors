from matplotlib import pyplot as plt
import numpy as np

nodeCount = []
calcTime = []

with open("time_for_octree.txt") as f:
   for line in f:
      values = line.split(" ")
      nodeCount.append(float(values[0]))
      calcTime.append(float(values[1]))

# plt.xscale("symlog", base=2)
plt.plot(nodeCount, calcTime, "-o")
plt.grid(True)
#plt.yticks(np.arange(0, 2.1, 0.1))
plt.xlabel("Количетсво точек квадратур")
plt.ylabel("Время построения, с")

plt.show()
