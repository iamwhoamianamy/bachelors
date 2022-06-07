from matplotlib import pyplot as plt
import numpy as np

quadTreeTime = []
pointTreeTime = []
pointOnQuadCount = []
m2mTime = []
m2ll2lTime = []

with open("fullFMM_by_quadratures.txt") as f:
   for line in f:
      values = line.split(" ")
      quadTreeTime.append(float(values[0]))
      pointTreeTime.append(float(values[1]))
      pointOnQuadCount.append(int(values[2]))
      m2mTime.append(float(values[3]))
      m2ll2lTime.append(float(values[4]))

totalPrep = [0] * len(pointOnQuadCount)
totalCalc = [0] * len(pointOnQuadCount)

for i in range(len(pointOnQuadCount)):
   totalPrep[i] = quadTreeTime[i] + m2mTime[i]
   totalCalc[i] = pointTreeTime[i] + m2ll2lTime[i]

plt.xscale("symlog", base=2)
plt.yscale("symlog", base=2)
# plt.plot(xs, withoutT)
# plt.plot(xs, withC)
# plt.plot(xs, withR)
plt.plot(pointOnQuadCount, totalPrep, "-o")
plt.plot(pointOnQuadCount, totalCalc, "-o")
plt.grid(True)
#plt.yticks(np.arange(0, 46, 2))
#plt.yticks(np.arange(0, 9, 1))
#plt.xlabel("Максимальное количество квадратур в листе")
plt.xlabel("Количество точек на количество квадратур")
plt.ylabel("Время выполнения, с")
plt.legend(["Построение первого дерева и M2M", 
            "Построение второго дерева, M2L и L2L"])

plt.show()
