%% Linear piston problem solver

clear all;
close all;

% Number of cells in flow domain
N = 64;

% Structural mass and stiffness
m = 2;
k = 1;

% Time step size and number of time steps
dt    = 0.1;
Ndt   = 100;

% Integration method:
%   theta = 0   : first order explicit Euler
%   theta = 1/2 : second order trapezoidal rule
%   theta = 1   : first order implicit Euler
theta = 0.0;

% spatial discretization:
% 
% dW/dt = [ As  Asf ] [Ws]
%         [ Afs Af  ] [Wf]
[As Asf Afs Af] = matrices(N,m,k);

% Initial condition: piston displaced unit 1
%                    based on exact solution
omega = exact_omega(m,k);
W0    = exact_sol(omega,N,0);

% Energy matrix: Energy = W' * E * W
E       = 0.5 * eye(2*N+2) / N;
E(1,1)  = 0.5 * m;
E(2,2)  = 0.5 * k;
E0      = W0' * E * W0;

% Monolithic
W_mono = W0;
A_mono = [ As Asf 
           Afs Af ];
L_mono = eye(2*N+2) - dt * theta     * A_mono;
R_mono = eye(2*N+2) + dt * (1-theta) * A_mono;
M_mono = inv(L_mono)*R_mono;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FILL IN THE PARTITIONG SCHEMES BELOW      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partitioned sequential Structure-Fluid
W_seqsf = W0;
L_seqsf = 
R_seqsf = 
M_seqsf = inv(L_seqsf)*R_seqsf;

% Partitioned sequential Fluid-Structure
W_seqfs = W0;
L_seqfs = 
R_seqfs = 
M_seqfs = inv(L_seqfs)*R_seqfs;

% Partitioned parallel
W_par = W0;
L_par = 
R_par = 
M_par = inv(L_par)*R_par;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tvec = (0:Ndt)*dt;
uvec = ones(5,1)*W0(1);
qvec = ones(5,1)*W0(2);
pvec = ones(5,1)*W0(N+2);
evec = zeros(5,1);

% Actual integration
for i=1:Ndt
    t       = i*dt;
    W_exact = exact_sol(omega,N,t);
    W_mono  = M_mono  * W_mono;
    W_seqsf = M_seqsf * W_seqsf;
    W_seqfs = M_seqfs * W_seqfs;
    W_par   = M_par   * W_par;
    
    uvec(:,i+1) = [W_exact(1)
                   W_mono(1) 
                   W_seqsf(1) 
                   W_seqfs(1) 
                   W_par(1)  ];
    qvec(:,i+1) = [W_exact(2)
                   W_mono(2) 
                   W_seqsf(2) 
                   W_seqfs(2) 
                   W_par(2)  ];
    pvec(:,i+1) = [W_exact(N+2)
                   W_mono(N+2) 
                   W_seqsf(N+2) 
                   W_seqfs(N+2) 
                   W_par(N+2)  ];
    evec(:,i+1) = [0
                   W_mono'  * E * W_mono  / E0 - 1
                   W_seqsf' * E * W_seqsf / E0 - 1
                   W_seqfs' * E * W_seqfs / E0 - 1
                   W_par'   * E * W_par   / E0 - 1 ];
end

showdata = [1 1 1 1 1];
for i=1:5
    qvec(i,:) = qvec(i,:)*showdata(i);
    uvec(i,:) = uvec(i,:)*showdata(i);
    pvec(i,:) = pvec(i,:)*showdata(i);
    evec(i,:) = evec(i,:)*showdata(i);
end

figure(1)
hold off
title('Piston displacement');
hold on
plot(tvec,qvec);
legend('Exact','Monolithic','Sequential S->F','Sequential F->S','Parallel','Location','Best');
xlabel('Time');
ylabel('Piston displacement');
figure(2)
hold off
title('Piston velocity');
hold on
plot(tvec,uvec);
legend('Exact','Monolithic','Sequential S->F','Sequential F->S','Parallel','Location','Best');
xlabel('Time');
ylabel('Piston velocity');
figure(3)
hold off
title('Pressure at piston');
hold on
plot(tvec,pvec);
legend('Exact','Monolithic','Sequential S->F','Sequential F->S','Parallel','Location','Best');
xlabel('Time');
ylabel('Interface pressure');
figure(4)
hold off
title('System energy change (E - E_{exact})/E_0');
hold on
plot(tvec,evec);
legend('Exact','Monolithic','Sequential S->F','Sequential F->S','Parallel','Location','Best');
xlabel('Time');
ylabel('System energy change');
