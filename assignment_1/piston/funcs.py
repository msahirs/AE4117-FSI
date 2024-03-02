import numpy as np


def get_sys_mat(N,m,k):
    As = np.zeros((2,2))
    Asf = np.zeros((2,2*N))
    Afs = np.zeros((2*N,2))
    Af = np.zeros((2*N,2*N))

    As[0,1] = -k/m;
    As[1,0] = 1;

    Asf[0,-1] = 1/m;

    Afs[-1,0]= -N;

    Af[:N-2, N+1:2*N-1] -= np.eye(N-1)
    Af[1:N-1, N:2*N-2] += np.eye(N-1)
    Af[N:2*N-2, 1:N-1] -= np.eye(N-1)
    Af[N+1:2*N-1  ,:N-2] += np.eye(N-1)

    Af[0,N] -= 1;
    Af[N,0] += 1;

    Af[N-1,2*N-1] += 1;
    Af[2*N-1,N-1] -= 1;

    Af *= (0.5 * N)

    return As, Asf, Afs, Af