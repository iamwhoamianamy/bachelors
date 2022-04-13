import numpy as np
import torch

# matrix = torch.tensor(np.zeros((3, 3), dtype=complex))
# vec = torch.tensor(np.zeros((3, 1), dtype=complex))

# matrix_dev = matrix.to("cuda")
# vec_dev = vec.to("cuda")

# torch.cuda.synchronize()

# res_dev = matrix_dev.matmul(vec)

# res = res_dev.to(torch.device("cpu"))

print(torch.cuda.is_available())