from funcs import get_sys_mat
from matplotlib import pyplot as plt
import numpy as np

N= 4

As, Asf, Afs, Af = get_sys_mat(N,1,1)

plt.title(f"n = {N}")
plt.imshow(Af)
plt.colorbar()

plt.show()