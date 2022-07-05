from matplotlib import pyplot as plt
import numpy as np

NM = []
timeInt = []
timeMatrices = []
timeFullFMM = []

with open("full_N_M_time.txt") as f:
   for line in f:
      values = line.split(" ")
      NM.append(float(values[0]))
      timeInt.append(float(values[1]))
      timeMatrices.append(float(values[2]))
      timeFullFMM.append(float(values[3]))

plt.xscale("symlog", base=2)
plt.yscale("symlog", base=2)
plt.plot(NM, timeInt, "-o")
plt.plot(NM, timeMatrices, "-o")
plt.plot(NM, timeFullFMM, "-o")
plt.grid(True)
plt.xlabel("Количество точек на количество точек квадратур")
plt.ylabel("Время построения, с")
plt.legend(["Интегрирование квадратурами", 
            "БММ с одним деревом",
            "БММ с двумя деревьями"])
plt.show()
