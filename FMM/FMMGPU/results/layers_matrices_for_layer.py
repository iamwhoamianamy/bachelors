from matplotlib import pyplot as plt
import numpy as np

xs = []
LKCPU = []
LACPU = []
MKCPU = []
MACPU = []
LKGPU = []
LAGPU = []
MKGPU = []
MAGPU = []

with open("layers_matrices_for_layer.txt") as f:
   for line in f.readlines():
      values = line.split(" ")
      LKCPU.append(float(values[0]))
      LACPU.append(float(values[1]))
      MKCPU.append(float(values[2]))
      MACPU.append(float(values[3]))
      LKGPU.append(float(values[4]))
      LAGPU.append(float(values[5]))
      MKGPU.append(float(values[6]))
      MAGPU.append(float(values[7]))
      xs.append(values[8])

start = 0
stop = 5

fig, axs = plt.subplots(1, 2)

#plt.yscale("symlog", base=10)
axs[0].plot(xs[start:stop], LKCPU[start:stop], "--b")
axs[0].plot(xs[start:stop], LACPU[start:stop], "-b")
axs[0].plot(xs[start:stop], MKCPU[start:stop], "--r")
axs[0].plot(xs[start:stop], MACPU[start:stop], "-r")
axs[0].plot(xs[start:stop], LKGPU[start:stop], "--c")
axs[0].plot(xs[start:stop], LAGPU[start:stop], "c")
axs[0].plot(xs[start:stop], MKGPU[start:stop], "--m")
axs[0].plot(xs[start:stop], MAGPU[start:stop], "m")
axs[0].grid(True)


start = stop - 1
stop = 12


axs[1].plot(xs[start:stop], LKCPU[start:stop], "--b")
axs[1].plot(xs[start:stop], LACPU[start:stop], "-b")
axs[1].plot(xs[start:stop], MKCPU[start:stop], "--r")
axs[1].plot(xs[start:stop], MACPU[start:stop], "-r")
axs[1].plot(xs[start:stop], LKGPU[start:stop], "--c")
axs[1].plot(xs[start:stop], LAGPU[start:stop], "c")
axs[1].plot(xs[start:stop], MKGPU[start:stop], "--m")
axs[1].plot(xs[start:stop], MAGPU[start:stop], "m")
axs[1].grid(True)


# axs[0].pl ("Узлов на слое")
axs[0].set_ylabel("Время выполнения, с")

for i in range(2):
   axs[i].legend(["LKCPU",
                  "LTCPU",
                  "MKCPU",
                  "MTCPU",
                  "LKGPU",
                  "LTGPU",
                  "MKGPU",
                  "MTGPU"])
   
   axs[i].set_xlabel("Количество узлов на слое")

plt.show()
