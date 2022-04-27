from matplotlib import pyplot as plt
import numpy as np

xs = [15728640]
integr = []
withoutT = []
withC = []
withR = []
withLCPU = []
withLGPU = []
withMCPU = []
withMGPU = []

with open("time for all.txt") as f:
   for line in f:
      values = line.split(" ")
      integr.append(float(values[1]))
      withoutT.append(float(values[2]))
      withC.append(float(values[3]))
      withR.append(float(values[4]))
      withLCPU.append(float(values[5]))
      withLGPU.append(float(values[6]))
      withMCPU.append(float(values[7]))
      withMGPU.append(float(values[8]))
      

for i in range(1, len(withC)):
   xs.append(xs[i - 1] * 2)
   
plt.xscale("symlog", base=2)
plt.yscale("symlog", base=2)
plt.plot(xs, withoutT, "-o")
plt.plot(xs, withC, "-o")
plt.plot(xs, withR, "-o")
plt.plot(xs, withLCPU, "-o")
plt.plot(xs, withLGPU, "-o")
plt.plot(xs, withMCPU, "-o")
plt.plot(xs, withMGPU, "-o")
plt.grid(True)
plt.xlabel("Количество точек на количество квадратур")
plt.ylabel("Время выполнения, с")
plt.legend(["БММ без переноса мультиполей", 
            "БММ с комплексным переносом", 
            "БММ с вещественным переносом", 
            "БММ с переносом по слоям, CPU", 
            "БММ с переносом по слоям, GPU", 
            "БММ с переносом матрицами, CPU", 
            "БММ с переносом матрицами, GPU"])
plt.show()
