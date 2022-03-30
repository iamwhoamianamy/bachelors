from matplotlib import pyplot as plt
import numpy as np

xs = [4]
withoutT = []
withT = []

with open("without with.txt") as f:
   for line in f:
      values = line.split(" ")
      withoutT.append(float(values[1]))
      withT.append(float(values[2]))

for i in range(1, len(withT)):
   xs.append(xs[i - 1] * 2)
   
plt.xscale("symlog", base=2)
plt.plot(xs, withoutT)
plt.plot(xs, withT)
plt.grid(True)
plt.yticks(np.arange(0, 10, 0.5))
plt.xlabel("Максимальное количество квадратур в листе")
plt.ylabel("Время выполнения, с")
plt.legend(["БММ без переноса мультиполей", "БММ с переносом мультиполей"])
plt.show()
