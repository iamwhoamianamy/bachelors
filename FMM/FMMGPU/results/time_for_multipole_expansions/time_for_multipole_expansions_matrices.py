from matplotlib import pyplot as plt
import numpy as np

xs = data = np.loadtxt("node_count.txt", dtype=np.float32, skiprows=0, usecols=(0))
CPU5 = data = np.loadtxt("time_for_multipole_expansions_order5_cpu.txt", dtype=np.float32, skiprows=0, usecols=(0))
GPU5 = data = np.loadtxt("time_for_multipole_expansions_order5_gpu.txt", dtype=np.float32, skiprows=0, usecols=(0))

CPU10 = data = np.loadtxt("time_for_multipole_expansions_order10_cpu.txt", dtype=np.float32, skiprows=0, usecols=(0))
GPU10 = data = np.loadtxt("time_for_multipole_expansions_order10_gpu.txt", dtype=np.float32, skiprows=0, usecols=(0))

CPU15 = data = np.loadtxt("time_for_multipole_expansions_order15_cpu.txt", dtype=np.float32, skiprows=0, usecols=(0))
GPU15 = data = np.loadtxt("time_for_multipole_expansions_order15_gpu.txt", dtype=np.float32, skiprows=0, usecols=(0))

n = len(xs)

CPU5 = CPU5[:n]
GPU5 = GPU5[:n]
CPU10 = CPU10[:n]
GPU10 = GPU10[:n]

#plt.xscale("symlog", base=2)
# plt.plot(xs, withoutT)
# plt.plot(xs, withC)
# plt.plot(xs, withR)
plt.plot(xs, CPU5, "--r")
plt.plot(xs, GPU5, "--m")
plt.plot(xs, CPU10, "--y")
plt.plot(xs, GPU10, "--b")
plt.plot(xs, CPU15, "--g")
plt.plot(xs, GPU15, "--c")
plt.grid(True)
#plt.yticks(np.arange(0, 46, 2))
#plt.yticks(np.arange(0, 9, 1))
#plt.xlabel("Максимальное количество квадратур в листе")
plt.xlabel("Количество узлов в дереве")
plt.ylabel("Время выполнения, с")

plt.legend(["CPU, максимальный порядок гармоник = 5", 
            "GPU, максимальный порядок гармоник = 5",
            "CPU, максимальный порядок гармоник = 10",
            "GPU, максимальный порядок гармоник = 10",
            "CPU, максимальный порядок гармоник = 15",
            "GPU, максимальный порядок гармоник = 15"])
plt.show()
