from matplotlib import pyplot as plt
import numpy as np

cubeHeight = []
m2mTime = []
m2ll2lTime = []

with open("time_by_cube_height.txt") as f:
   for line in f:
      values = line.split(" ")
      cubeHeight.append(int(values[2]) * 2)
      m2mTime.append(float(values[4]))
      m2ll2lTime.append(float(values[5]))

#plt.xscale("symlog", base=2)
#plt.yscale("symlog", base=2)
# plt.plot(xs, withoutT)
# plt.plot(xs, withC)
# plt.plot(xs, withR)
plt.plot(cubeHeight, m2mTime, "-o")
plt.plot(cubeHeight, m2ll2lTime, "-o")
plt.grid(True)
#plt.yticks(np.arange(0, 46, 2))
#plt.yticks(np.arange(0, 9, 1))
#plt.xlabel("Максимальное количество квадратур в листе")
plt.xlabel("Высота параллелепипеда, м")
plt.ylabel("Время выполнения, с")
plt.legend(["Построение первого дерева и M2M", 
            "Построение второго дерева, M2L и L2L"])

plt.show()
