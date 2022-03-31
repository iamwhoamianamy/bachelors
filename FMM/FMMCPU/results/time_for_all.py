from matplotlib import pyplot as plt
import numpy as np

xs = [15728640]
integr = []
withoutT = []
withT = []

with open("time for all.txt") as f:
   for line in f:
      values = line.split(" ")
      integr.append(float(values[1]))
      withoutT.append(float(values[2]))
      withT.append(float(values[3]))

for i in range(1, len(withT)):
   xs.append(xs[i - 1] * 2)
   
plt.xscale("symlog", base=2)
plt.yscale("symlog", base=2)
plt.plot(xs, integr, "-o")
plt.plot(xs, withoutT, "-o")
plt.plot(xs, withT, "-o")
plt.grid(True)
plt.xlabel("Количество точек на количество квадратур")
plt.ylabel("Время выполнения, с")
plt.legend(["Интегрирование квадратурами", "БММ без переноса мультиполей", "БММ с переносом мультиполей"])
plt.show()
