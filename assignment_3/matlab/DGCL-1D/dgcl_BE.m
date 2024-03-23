close all
clear all
%% Input data
N  = 20;        % number of cells

dt   = 0.1;     % time step
tend = 10;      % total simulation time

fmt = '-b';     % linecolor for plot of final solution

show_sol = 1;   % (=1) show intermediate solutions during simulation

method = 3;     % method = 1 : exact at t(n+1)
                % method = 2 : exact at t(n+1/2)
                % method = 3 : DGCL

dt_values = [0.1, 0.05, 0.025];
for dt = dt_values
  if dt==0.1
    fmt = '-b';
  elseif dt==0.05
    fmt = '-r';
  else
    fmt = '-g';
  end

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
  dxidt_dgcl    = (xi - xi_tn) / dt;                            % face velocity satisfying D-GCL
  
  if (method == 1)
    dxidt = dxidt_exnp1;
  elseif (method ==2)
    dxidt = dxidt_exnp1_2;
  else
    dxidt = dxidt_dgcl;
  end
    
  % Setting up system to solve; for each cell:
  %
  %   (u*dx)^(n+1) - (u*dx)^n      
  %   ------------------------- - [ u_leftface * dxidt_leftface * (-1) + u_rightface * dxidt_rightface * (+1) ]^(n+1) = 0
  %             dt                 
  %
  % which can be written for the total system as:
  %
  % L u^(n+1) = (u*dx)^n
  %
  % We use a simple avarage to compute the solution at a cell face:
  %   u_face = 0.5 * (u_leftcell + u_rightcell)
  %
  % Internal cells
  for i=2:N-1
    L(i,i-1:i+1) = [0 dx(i) 0] - dt * 0.5*[-dxidt(i) dxidt(i+1)-dxidt(i) dxidt(i+1)];
  end
  % Boundary cells
  L(1,1:2)   = [dx(1) 0] - dt * 0.5*[dxidt(2)-2*dxidt(1) dxidt(2)];
  L(N,N-1:N) = [0 dx(N)] - dt * 0.5*[-dxidt(N) 2*dxidt(N+1)-dxidt(N)];

  % solve system
  u = (L\(dx_tn.*u_tn)')';

  % Show intermediate solution
  if (show_sol)
    figure(1);
    hold off;
    plot(x,u,'-x');
    axis([0 1 0.8 1.2]);
%    pause
  end
  
end

%% Plot final solution
figure(2);
hold on;
title(['Solution at t=' num2str(tend)]);
ylabel('Solution');
xlabel('x');
plot(x,u,fmt);
legend(['dt=0.1'], ['dt=0.05'], ['dt=0.025']);


end


saveas(figure(2), 'method_3_q21.png');
hold off;
