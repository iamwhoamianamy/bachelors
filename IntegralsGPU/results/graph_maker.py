import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def read_values(device, datatype):
   values = []
   with open("laplace_solver_structs_results_" + device + "_" + datatype + ".txt") as f:
      for line in f:
         values.append(float(line))
   return values

def fix_cpu_values_double(values):
   for i in range(len(values)):
      values[i] /= 32
   return values

def fix_cpu_values_float(values):
   for i in range(len(values)):
      values[i] /= 64
   return values

def plot_data(xs, values, device, datatype):
   plt.xlabel("Количество точек")
   plt.ylabel("Время выполнения, c")
   plt.xticks(np.arange(0, xs[-1] * 2, xs[-1] / 8))
   plt.yticks(np.arange(0, values[-1] * 2, values[-1] / 16))
   #plt.plot(xs[:n_CPU], ys_CPU, '-ro')
   
   color = 'b' if device == "GPU" else 'r'
   stroke = '--' if datatype == "double" else '-'
   
   plt.plot(xs, values, stroke + color + 'o')
   plt.grid(True)
   
def plot_one(device, datatype):
   values = []

   values = read_values(device, datatype)

   if device == "CPU":
      if datatype == "double":
         values = fix_cpu_values_double(values)
      else:
         values = fix_cpu_values_float(values)
         

   count = len(values)

   xs = []

   for i in range(1, count + 1):
      xs.append(2**i)

   plot_data(xs, values, device, datatype)
   
# plot_one("GPU", "float")
# plot_one("GPU", "double")

# plt.legend(["CPU float", "CPU double"])

plot_one("CPU", "float")
plot_one("CPU", "double")
plt.legend(["CPU float", "CPU double"])

plt.show()
