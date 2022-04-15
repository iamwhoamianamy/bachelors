import timeit
import numpy as np
import torch
import cupy as cnp

order = 10
harmonicLength = (order + 1)**2

def matrix_from_file(filename):
   res = []
   
   with open(filename) as f:
      y = 0
      for line in f.readlines():
         temp = []
         
         for token in line.split(" "):
            if token != "\n":
               token_splited = token.replace("(", "").replace(")", "").split(',')
               temp.append(complex(np.single(token_splited[0]), np.single(token_splited[1])))
               
         res.append(temp)
            
   return res
   
def calc_for_one_layer_one_orientation(regular, expansions):
   regular_gpu = regular.to("cuda")
   expansion_x_gpu = expansions[0].to("cuda")
   expansion_y_gpu = expansions[1].to("cuda")
   expansion_z_gpu = expansions[2].to("cuda")
   torch.cuda.synchronize()
   
   start = timeit.default_timer()
   
   x_contribution_gpu = regular_gpu.mm(expansion_x_gpu)
   y_contribution_gpu = regular_gpu.mm(expansion_y_gpu)
   z_contribution_gpu = regular_gpu.mm(expansion_z_gpu)
   torch.cuda.synchronize()
   
   stop = timeit.default_timer()
   
   x_contribution = x_contribution_gpu.to("cpu")
   y_contribution = y_contribution_gpu.to("cpu")
   z_contribution = z_contribution_gpu.to("cpu")

   print((stop - start) * 8)

   return [x_contribution, y_contribution, z_contribution]

def calc_for_one_layer_one_orientation_cupy(regular, expansions):
   regular_gpu = cnp.array(regular, dtype=np.complex64)
   expansion_x_gpu = cnp.array(expansions[0], dtype=np.complex64)
   expansion_y_gpu = cnp.array(expansions[1], dtype=np.complex64)
   expansion_z_gpu = cnp.array(expansions[2], dtype=np.complex64)
      
   start = timeit.default_timer()
   
   contribution_x_gpu = cnp.matmul(regular_gpu, expansion_x_gpu)
   contribution_y_gpu = cnp.matmul(regular_gpu, expansion_y_gpu)
   contribution_z_gpu = cnp.matmul(regular_gpu, expansion_z_gpu)

   stop = timeit.default_timer()

   print((stop - start) * 8)
   
#regular_0 = torch.tril(matrix_from_file("regular_0.txt"))
# expansion_x = torch.tensor(matrix_from_file("expansion_x.txt"), dtype=torch.complex64) 
# expansion_y = torch.tensor(matrix_from_file("expansion_y.txt"), dtype=torch.complex64) 
# expansion_z = torch.tensor(matrix_from_file("expansion_z.txt"), dtype=torch.complex64) 

regular_0 = matrix_from_file("regular_0.txt")
expansion_x = matrix_from_file("expansion_x.txt")
expansion_y = matrix_from_file("expansion_y.txt")
expansion_z = matrix_from_file("expansion_z.txt")

print("done reading!")
print("expansions count: " + str(len(expansion_x[0])))

calc_for_one_layer_one_orientation_cupy(regular_0, [expansion_x, expansion_y, expansion_z])