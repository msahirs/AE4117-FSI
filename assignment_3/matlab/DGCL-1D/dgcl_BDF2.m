%% Input data
N  = 20;        % number of cells

dt   = 0.1;     % time step
tend = 100;      % total simulation time

fmt = 'b-';     % linecolor for plot of final solution

show_sol = 1;   % (=1) show intermediate solutions during simulation

method = 3;     % method = 1 : exact at t(n+1)
                % method = 2 : exact at t(n+1/2)
                % method = 3 : DGCL
if (method == 1)
  disp('Mesh velocities are EXACT at t(n+1)');
elseif (method == 2)
  disp('Mesh velocities are EXACT at t(n+1/2)');
else
  disp('Mesh velocities satisfy DGCL');
end

%% Initialisation

dx = zeros(1,N);    % cell volumes
xi = zeros(1,N+1);  % face centers
x  = zeros(1,N);    % cell centers
u  = ones(1,N);     % solution

dx_tn = zeros(1,N);    % cell volumes at tn
xi_tn = zeros(1,N+1);  % face centers at tn
u_tn  = ones(1,N);     % solution at tn

dx_tnm1 = zeros(1,N);    % cell volumes at tn-1
xi_tnm1 = zeros(1,N+1);  % face centers at tn-1
u_tnm1  = ones(1,N);     % solution at tn-1

xi  = (0:N)/N;              % face centers
xi0 = (0:N)/N;              % face centers at t=0
x   = (1:N)/N - 0.5/N;      % cell centers
dx  = xi(2:N+1) - xi(1:N);  % cell volumes

% Plot initial solution
figure(1);
hold off;
axis([0 1 0.98 1.02]);
plot(x,u,'-x');
hold off;

L = zeros(N);   % Discretization matrix

alpha_BE = [1 -1 0]; % Time discretization coefficients (first step Backward Euler)
alpha_BDF = [1.5 -2 0.5]; % Time discretization coefficients (first step Backward Euler)
alpha = alpha_BE; %[1 -1 0]; % Time discretization coefficients (first step Backward Euler)

%% Simulation loop
for t=dt:dt:tend

  dx_tn = dx;  % store cell volume at tn
  xi_tn = xi;  % store face center at tn
  u_tn  = u;   % store solution at tn
  
  xi     = xi0 + sin(2*pi*t) * sin(2*pi*xi0) / N;   % face centers at tn+1
  dx     = xi(2:N+1)-xi(1:N);                       % cell volumes at tn+1
  x      = xi(1:N) + dx/2;                          % cell centers at tn+1
  
  dxidt_exnp1   = 2*pi*cos(2*pi*(t))*sin(2*pi*xi0) / N;         % exact face velocity at tn+1
  dxidt_exnp1_2 = 2*pi*cos(2*pi*(t-dt/2))*sin(2*pi*xi0) / N;    % exact face velocity at tn+1/2

  %% IMPLEMENTATION OF DGCL
  dxidt_dgcl    = (alpha(1) * xi + alpha(2) * xi_tn + alpha(3) * xi_tnm1) / dt; % face velocity satisfying D-GCL
  %%

  if (method == 1)
    dxidt = dxidt_exnp1;
  elseif (method ==2)
    dxidt = dxidt_exnp1_2;
  else
    dxidt = dxidt_dgcl;
  end
    
  % Setting up system to solve; for each cell:
  %
  %    i=3            (u*dx)^(n+2-i)      
  %   SUM   alpha(i) -------------- - [ u_leftface * dxidt_leftface * (-1) + u_rightface * dxidt_rightface * (+1) ] = 0
  %    i=1                 dt                 
  %
  % which can be written for the total system as:
  %
  % L u^(n+1) = alpha(0) * (u*dx)^n + alpha(-1) * (u*dx)^(n-1)
  %
  % We use a simple avarage to compute the solution at a cell face:
  %   u_face = 0.5 * (u_leftcell + u_rightcell)
  %
  % Internal cells
  for i=2:N-1
    L(i,i-1:i+1) = [0 alpha(1)*dx(i) 0] - dt * 0.5*[-dxidt(i) dxidt(i+1)-dxidt(i) dxidt(i+1)];
  end
  % Boundary cells
  L(1,1:2)   = [alpha(1)*dx(1) 0] - dt * 0.5*[dxidt(2)-2*dxidt(1) dxidt(2)];
  L(N,N-1:N) = [0 alpha(1)*dx(N)] - dt * 0.5*[-dxidt(N) 2*dxidt(N+1)-dxidt(N)];

  % solve system
  u = (L\(-alpha(2)*(dx_tn.*u_tn)'-alpha(3)*(dx_tnm1.*u_tnm1)'))';

  % update old variables
  dx_tnm1 = dx_tn;
  u_tnm1  =  u_tn;
  xi_tnm1 = xi_tn;
  
  % After first time step we can use 
  alpha= alpha_BDF;
  
  % Show intermediate solution
%   if (show_sol)
%     figure(1);
%     hold off;
%     plot(x,u,'-x');
%     axis([0 1 0.8 1.2]);
% %    pause
%   end
  
end

%% Plot final solution
figure(2);
hold on;
title(['Solution at t=' num2str(tend)]);
ylabel('Solution');
xlabel('x');
plot(x,u,fmt);
