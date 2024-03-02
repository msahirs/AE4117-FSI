import numpy as np


def get_sys_mat(N,m,k):
    As = np.zeros((2,2))
    Asf = np.zeros((2,2*N))
    Afs = np.zeros((2*N,2))
    Af = np.zeros((2*N,2*N))

    As[0,1] = -k/m;
    As[1,0] = 1;

    Asf[0,-N-1] = 1/m;

    Afs[-1,0]= -N;

    Af[:N-1, N+1:2*N] -= np.eye(N-1)
    Af[1:N, N:2*N-1] += np.eye(N-1)
    Af[N:2*N-1, 1:N] -= np.eye(N-1)
    Af[N+1:2*N  ,:N-1] += np.eye(N-1)

    Af[0,N] -= 1;
    Af[N,0] += 1;

    Af[N-1,2*N-1] += 1;
    Af[2*N-1,N-1] -= 1;

    Af *= (0.5 * N)

    return As, Asf, Afs, Af

def get_exact_initsol(o,N,t):
    W_s = np.array([-o*np.sin(o*t), np.cos(o*t) ])
    W_f = np.zeros((2*N,1));
    x = np.arange(0.5/N,1/N,(1-0.5/N))
    W_f[:N] = - o/np.sin(o) * np.cos(o*x[:,np.newaxis]) * np.cos(o*t)
    W_f[N:2*N] = - o/np.sin(o) * np.sin(o*x[:,np.newaxis]) * np.sin(o*t)

    return W_f, W_s

def get_exact_sol(o,N,t):
    x = np.arange(0.5/N,1/N,(1-0.5/N))
    W = np.zeros((2*N+2,len(t)));
    W[:2,:]       = np.array([ -o*np.sin(o*t),
                        np.cos(o*t) ])
    W[2:N+2,:]= - o / np.sin(o) * np.cos(o*x[:,np.newaxis]) * np.cos(o*t)
    W[N+1:2*N+2,:] = - o / np.sin(o) * np.sin(o*x[:,np.newaxis]) * np.sin(o*t)

    return W


def get_exact_omega(m,k):
    nmax_iter = 100
    phase_min = 1e-10
    phase_max = np.pi
    npi       = 0
    phase     = 0.5 * np.pi
    o_old     = -1
    o         =  1

    i = 0

    while not((o == o_old) or (i > nmax_iter)):
        o_old = o
        i+=1
        o = npi * np.pi + phase_min
        err_min = ((m * o*o - k) * np.sin(phase_min) - o * np.cos(phase_min)) \
                (k + o + m * o*o)
        err_max = ((m * o*o - k) * np.sin(phase_max) - o * np.cos(phase_max)) \
                (k + o + m * o*o)
        o = npi * np.pi + phase
        err = ((m * o*o - k) * np.sin(phase) - o * np.cos(phase)) \
            (k + o + m * o*o)
        
        if (err_min * err_max > 0):
            raise BaseException("ERRRORRRRR")
        elif (err_min * err_max < 0):
            phase_max = phase
            phase = phase_min + (phase - phase_min) / (err_min - err) * err_min
        else:
            phase_min = phase;
            phase = phase + (phase_max - phase) / (err - err_max) * err

    if ( o_old != o ):
        print('Could not find correct mode...')

    return o