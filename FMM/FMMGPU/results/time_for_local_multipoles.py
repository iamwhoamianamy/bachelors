from matplotlib import pyplot as plt
import numpy as np

xs = [8]
withoutT = []
withC = []
withR = []
withLCPU = []
withLGPU = []
withMCPU = []
withMGPU = []


with open("time_for_local_multipoles.txt") as f:
   for line in f:
      values = line.split(" ")
      withoutT.append(float(values[1]))
      withC.append(float(values[2]))
      withR.append(float(values[3]))
      withLCPU.append(float(values[4]))
      withLGPU.append(float(values[5]))
      withMCPU.append(float(values[6]))
      withMGPU.append(float(values[7]))

for i in range(1, len(withC)):
   xs.append(xs[i - 1] * 2)
   
plt.xscale("symlog", base=2)
plt.plot(xs, withoutT)
plt.plot(xs, withC)
plt.plot(xs, withR)
plt.plot(xs, withLCPU, "--b")
plt.plot(xs, withLGPU, "--c")
plt.plot(xs, withMCPU, "--r")
plt.plot(xs, withMGPU, "--m")
plt.grid(True)
plt.yticks(np.arange(0, 45, 2))
#plt.yticks(np.arange(0, 9, 1))
plt.xlabel("Максимальное количество квадратур в листе")
plt.ylabel("Время выполнения, с")
plt.legend(["Без переноса мультиполей", 
            "Комплексный перенос", 
            "Вещественный перенос", 
            "Перенос по слоям, CPU", 
            "Перенос по слоям, GPU", 
            "Перенос матрицами, CPU", 
            "Перенос матрицами, GPU"])

# plt.legend(["Перенос по слоям, CPU", 
#             "Перенос по слоям, GPU", 
#             "Перенос матрицами, CPU", 
#             "Перенос матрицами, GPU"])
plt.show()
