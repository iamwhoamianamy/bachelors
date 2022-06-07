from matplotlib import pyplot as plt
import numpy as np

xs = []
withoutT = []
withC = []
withR = []
withLCPU = []
withLGPU = []
withMCPU = []
withMGPU = []

with open("time_for_local_multipoles_nodes.txt") as f:
   for line in f:
      values = line.split(" ")
      xs.append(float(values[0]) / 2)
      withoutT.append(float(values[1]))
      withC.append(float(values[2]))
      withR.append(float(values[3]))
      withLCPU.append(float(values[4]))
      withLGPU.append(float(values[5]))
      withMCPU.append(float(values[6]))
      withMGPU.append(float(values[7]))

#plt.xscale("symlog", base=2)
# plt.plot(xs, withoutT)
# plt.plot(xs, withC)
# plt.plot(xs, withR)
plt.plot(xs, withLCPU, "--b")
plt.plot(xs, withLGPU, "--c")
plt.plot(xs, withMCPU, "--r")
plt.plot(xs, withMGPU, "--m")
plt.grid(True)
plt.yticks(np.arange(0, 46, 2))
#plt.yticks(np.arange(0, 9, 1))
#plt.xlabel("Максимальное количество квадратур в листе")
plt.xlabel("Количество узлов в дереве")
plt.ylabel("Время выполнения, с")
# plt.legend(["Без переноса мультиполей", 
#             "Комплексный перенос", 
#             "Вещественный перенос", 
#             "Перенос по слоям, CPU", 
#             "Перенос по слоям, GPU", 
#             "Перенос матрицами, CPU", 
#             "Перенос матрицами, GPU"])

plt.legend(["Перенос по слоям, CPU", 
            "Перенос по слоям, GPU", 
            "Перенос матрицами, CPU", 
            "Перенос матрицами, GPU"])
plt.show()
