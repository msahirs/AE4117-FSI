import numpy as np
import funcs
from matplotlib import pyplot as plt

# Number of cells in flow domain
N = 64

# Structural mass and stiffness
m = 2
k = 1


# Time step size and number of time steps
dt    = 0.1
Ndt   = 1000

# Integration method:
#   theta = 0   : first order explicit Euler
#   theta = 1/2 : second order trapezoidal rule
#   theta = 1   : first order implicit Euler
theta = 0.5

# spatial discretization:
# 
# dW/dt = [ As  Asf ] [Ws]
#         [ Afs Af  ] [Wf]
As, Asf, Afs, Af = funcs.get_sys_mat(N,m,k)

# Initial condition: piston displaced unit 1
#                    based on exact solution
omega = funcs.get_exact_omega(m,k)

W0    = funcs.get_exact_sol(omega,N,0)
# print(W0)


# Energy matrix: Energy = W' * E * W
E       = 0.5 * np.eye(2*N+2) / N
E[0,0]  = 0.5 * m
E[1,1]  = 0.5 * k
E0      = W0.T @ E @ W0


# Monolithic
W_mono = np.copy(W0)
A_mono = np.vstack((np.hstack((As, Asf)),
           np.hstack((Afs, Af))))
L_mono = np.eye(2*N+2) - dt * theta     * A_mono
R_mono = np.eye(2*N+2) + dt * (1-theta) * A_mono
M_mono = np.linalg.inv(L_mono)@R_mono

###############################################
## FILL IN THE PARTITIONG SCHEMES BELOW      ##
###############################################
# print(np.block([[As, 0], [Afs, Af]]))
# Partitioned sequential Structure-Fluid
W_seqsf = np.copy(W0)
L_seqsf = np.eye(2*N+2) - dt * theta* np.block([[As, np.zeros((2,2*N))],
                                                    [Afs, Af]])
R_seqsf = np.eye(2*N+2) + dt * (1-theta)     * np.block([[np.zeros((2,2)), Asf], 
                                                    [np.zeros((2*N,2)), np.zeros((2*N,2*N))]])
M_seqsf = np.linalg.inv(L_seqsf) @ R_seqsf

# Partitioned sequential Fluid-Structure
W_seqfs = np.copy(W0)
L_seqfs = np.eye(2*N+2) - dt * theta     * np.block([[As, Asf], 
                                                    [np.zeros((2*N,2)), Af]])
R_seqfs = np.eye(2*N+2) + dt * (1-theta)     * np.block([[np.zeros((2,2)), np.zeros((2,2*N))], 
                                                    [Afs, np.zeros((2*N,2*N))]])
M_seqfs = np.linalg.inv(L_seqfs) @ R_seqfs

# Partitioned parallel
W_par = np.copy(W0)
L_par = np.eye(2*N+2) - dt * theta     * np.block([[As, np.zeros((2,2*N))], 
                                                    [np.zeros((2*N,2)), Af]])
R_par = np.eye(2*N+2) + dt * (1-theta)     * np.block([[np.zeros((2,2)), Asf], 
                                                    [Afs, np.zeros((2*N,2*N))]])
M_par = np.linalg.inv(L_par) @ R_par
###############################################

tvec = np.arange(Ndt)*dt
uvec = np.ones((5,Ndt))
qvec = np.ones((5,Ndt))
pvec = np.ones((5,Ndt))
evec = np.zeros((5,Ndt))

# print(funcs.get_exact_sol(omega,N,0))


# Actual integration
for i in range(Ndt):

    t       = i*dt
    W_exact = funcs.get_exact_sol(omega,N,t)

    W_mono  = M_mono  @ W_mono
    W_seqsf = M_seqsf @ W_seqsf
    W_seqfs = M_seqfs @ W_seqfs
    W_par   = M_par   @ W_par

    uvec[:,i] = np.array([W_exact[0,0],
                   W_mono[0,0] ,
                   W_seqsf[0,0] ,
                   W_seqfs[0,0] ,
                   W_par[0,0]  ])
    qvec[:,i] = np.array([W_exact[1,0],
                   W_mono[1,0] ,
                   W_seqsf[1,0] ,
                   W_seqfs[1,0] ,
                   W_par[1,0]  ])
    
    pvec[:,i] = np.array([W_exact[N+1,0],
                   W_mono[N+1,0] ,
                   W_seqsf[N+1,0] ,
                   W_seqfs[N+1,0] ,
                   W_par[N+1,0]  ])
    
    
    evec[:,i] = np.array([0,
                   (W_mono.T  @ E @ W_mono  / E0 -1 )[0,0],
                   (W_seqsf.T @ E @ W_seqsf / E0 -1)[0,0],
                   (W_seqfs.T @ E @ W_seqfs / E0 -1)[0,0],
                   (W_par.T   @ E @ W_par   / E0 -1)[0,0] ])

# for i in range(5):
#     qvec[i,:] = qvec[i,:]*showdata[i]
#     uvec[i,:] = uvec[i,:]*showdata[i]
#     pvec[i,:] = pvec[i,:]*showdata[i]
#     evec[i,:] = evec[i,:]*showdata[i]

lst = [uvec, qvec, pvec, evec]
quote = ['u', 'q', 'p', 'e']
figure_name = ['theta_uvec', 'theta_qvec', 'theta_pvec', 'theta_evec']
for i in range(0,4):
    plt.figure(figure_name[i])
    plt.clf()
    # # funky truncation fix to preserve graph for forward euler, monolithic solution
    # if (theta == 0):
    #     truc_limit = 5
    #     plt.plot(tvec[:truc_limit],lst[i][1][:truc_limit], label = 'mono (divergent)') # plot only the first k few time-steps
    # else:
    #     plt.plot(tvec,lst[i][1],'--', dashes=(5, 3), label = 'mono') # plot entire time series
    plt.plot(tvec,lst[i][0], label = 'exact')
    plt.plot(tvec,lst[i][1],'--', dashes=(5, 3), label = 'mono')
    plt.plot(tvec,lst[i][2],  label = 's->f')
    plt.plot(tvec,lst[i][3],'--',dashes=(5, 1), label = 'f->s')
    plt.plot(tvec,lst[i][4], '--', dashes=(5, 5), label = 'parallel')
    plt.xlabel("Time [s]")
    plt.ylabel(f"{quote[i]}")
    plt.legend()
    plt.show()
# print(lst[3][1][-1]-lst[3][0][-1])
# print(lst[3][2][-1]-lst[3][0][-1])
# print(lst[3][3][-1]-lst[3][0][-1])
# print(lst[3][4][-1]-lst[3][0][-1])
print(f'Remaining energy in system: Partioned Scheme vs Monolithic \nSerial S->F {np.abs(evec[2, -1]- evec[1,-1])} \nSerial F->S {np.abs(evec[3, -1]- evec[1,-1])} \nParallel: {np.abs(evec[4, -1]- evec[1,-1])} ')
print(f'Remaining energy in system: Partioned Scheme vs Exact \nSerial S->F {np.abs(evec[2, -1]- evec[0,-1])} \nSerial F->S {np.abs(evec[3, -1]- evec[0,-1])} \nParallel: {np.abs(evec[4, -1]- evec[0,-1])} ')