from funcs import get_sys_mat, get_exact_sol
from matplotlib import pyplot as plt
import numpy as np

N= 4

print(get_exact_sol(1,N,np.linspace(0,1,1000)))

As, Asf, Afs, Af = get_sys_mat(N,1,1)



A = np.vstack((np.hstack((As,Asf)), np.hstack((Afs,Af))))

plt.title(f"n = {N}")
plt.imshow(A)
plt.colorbar()

plt.show()