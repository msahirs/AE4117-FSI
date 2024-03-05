import numpy as np
import funcs
from matplotlib import pyplot as plt
from scipy.signal import find_peaks


# Number of cells in flow domain
N = 64

# Structural mass and stiffness
m = 2
k = 1


# Time step size and number of time steps
dt    = 0.1
Ndt   = 100

# Integration method:
#   theta = 0   : first order explicit Euler
#   theta = 1/2 : second order trapezoidal rule
#   theta = 1   : first order implicit Euler
theta = 1

# spatial discretization:
# 
# dW/dt = [ As  Asf ] [Ws]
#         [ Afs Af  ] [Wf]
As, Asf, Afs, Af = funcs.get_sys_mat(N,m,k)

# Initial condition: piston displaced unit 1
#                    based on exact solution
omega = funcs.get_exact_omega(m,k)
# print(omega)
W0    = funcs.get_exact_sol(omega,N,0)


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
M_mono = np.linalg.inv(L_mono)*R_mono

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
uvec = np.ones((5,Ndt))*W0[0]
qvec = np.ones((5,Ndt))*W0[1]
pvec = np.ones((5,Ndt))*W0[N+1]
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
                   (W_mono.T  * E * W_mono  / E0 )[0,0],
                   (W_seqsf.T * E * W_seqsf.T / E0)[0,0],
                   (W_seqfs.T * E * W_seqfs.T / E0 )[0,0],
                   (W_par.T   * E * W_par.T   / E0)[0,0] ])

## if FE:
# coefficients = np.polyfit(tvec, evec[2], 2) # FE, S->F SERIAL
# p = np.poly1d(coefficients)



# for i in range(5):
#     qvec[i,:] = qvec[i,:]*showdata[i]
#     uvec[i,:] = uvec[i,:]*showdata[i]
#     pvec[i,:] = pvec[i,:]*showdata[i]
#     evec[i,:] = evec[i,:]*showdata[i]

plt.title(f"u v time, dt = {dt}, N = {N}")
plt.clf()
plt.plot(tvec,evec[0],'--',label = 'exact')
plt.plot(tvec,evec[1], label = 'mono')
plt.plot(tvec,evec[2], label = 's->f')
plt.plot(tvec,evec[3], label = 'f->s')
plt.plot(tvec,evec[4], label = 'parallel')
# plt.plot(tvec, p(tvec), ':', label = 'line fit for forward euler s->f = paralel')


plt.xlabel("Time [s]")
plt.ylabel("u")
plt.legend()

plt.show()
# figure[0]
# hold off
# title('Piston displacement')
# hold on
# plot(tvec,qvec)
# legend('Exact','Monolithic','Sequential S->F','Sequential F->S','Parallel','Location','Best')
# xlabel('Time')
# ylabel('Piston displacement')
# figure[1]
# hold off
# title('Piston velocity')
# hold on
# plot(tvec,uvec)
# legend('Exact','Monolithic','Sequential S->F','Sequential F->S','Parallel','Location','Best')
# xlabel('Time')
# ylabel('Piston velocity')
# figure(3)
# hold off
# title('Pressure at piston')
# hold on
# plot(tvec,pvec)
# legend('Exact','Monolithic','Sequential S->F','Sequential F->S','Parallel','Location','Best')
# xlabel('Time')
# ylabel('Interface pressure')
# figure(4)
# hold off
# title('System energy change (E - E_{exact})/E_0')
# hold on
# plot(tvec,evec)
# legend('Exact','Monolithic','Sequential S->F','Sequential F->S','Parallel','Location','Best')
# xlabel('Time')
# ylabel('System energy change')