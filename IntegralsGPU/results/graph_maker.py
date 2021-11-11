import matplotlib.pyplot as plt
import numpy as np

points_count = []
cpu_with_vector3s_t = []
gpu_with_vector3s_t = []

with open("laplace_solver_vector3s_test_results.txt") as f:
   for line in f:
      cols = line.split("\t")
      points_count.append(float(cols[0]))
      cpu_with_vector3s_t.append(float(cols[1]))
      gpu_with_vector3s_t.append(float(cols[2]))
   
n = len(cpu_with_vector3s_t)
x = np.arange(2, 2**n)

plt.figure(figsize=(20, 10))
plt.subplot(1, 2, 1)

plt.xticks(x)
# plt.xlim([cpu_with_vector3s_t[0], cpu_with_vector3s_t[n - 1]])
plt.xlabel("Количество точек")
plt.ylabel("Время выполнения, c")
plt.plot(points_count, cpu_with_vector3s_t)
plt.plot(points_count, gpu_with_vector3s_t)

plt.subplot(1, 2, 2)
plt.xlabel("Количество точек")
plt.ylabel("Ускорение видеокарты, c")
plt.plot(points_count, [float(cpu_with_vector3s_t[i])  / float(gpu_with_vector3s_t[i]) for i in range(len(cpu_with_vector3s_t))])
plt.show()