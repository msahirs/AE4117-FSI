close all
clear all
%% Input data
N  = 20;        % number of cells

dt   = 0.1;     % time step
tend = 10;      % total simulation time

fmt = '--b';     % linecolor for plot of final solution

show_sol = 0;   % (=1) show intermediate solutions during simulation
 
u_ex = ones(N); % unsure: creates exact scalar field solution

%% Run simulation for decreasing dt
k_max = 6;
e1_arr = zeros(k_max, 3);
e2_arr = zeros(k_max, 3);
dt_arr = zeros(k_max, 1); 
for method = 1:1:3
    fprintf('\nRunning Method %d', method)
    for k=1:1:k_max+1
         dt = 0.1*2^(-(k-1)); % time-step
         dt_arr(k, 1) = dt; 
    
         e1 = 0; % error 1: RMS error at final time-step
         e2 = 0; % error 2: RMS error at all time-steps
    
        %% Initialisation
        dx = zeros(1,N);    % cell volumes
        xi = zeros(1,N+1);  % face centers
        x  = zeros(1,N);    % cell centers
        u  = ones(1,N);     % solution
        
        xi  = (0:N)/N;              % face centers
        xi0 = (0:N)/N;              % face centers at t=0
        x   = (1:N)/N - 0.5/N;      % cell centers
        dx  = xi(2:N+1) - xi(1:N);  % cell volumes

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
        
          % error 1 computation
          if abs(t-tend) < 0.01*dt
              for i = 1:N
                e1 = e1 + (u(i)-u_ex(i))*(u(i)-u_ex(i));
              end
          end
        
          % error 2 computation
          for i=1:N
            e2 = e2 + (u(i)-u_ex(i))*(u(i)-u_ex(i));
          end

        end
        %% Error Computation
        e1 = sqrt(e1)/N;
        M = round(tend/dt);
        e2 = sqrt(e2)/(N*M);
    
        e1_arr(k, method) = e1;
        e2_arr(k, method) = e2;
    end
end

% error 1 plots 
figure(1);
hold("off");
loglog(dt_arr, e1_arr(:,1), 'b-x'); % method 1
hold("on");
loglog(dt_arr, e1_arr(:,2), 'r-x'); % method 2
loglog(dt_arr, e1_arr(:,3), 'g-x'); % method 3
legend('Method 1: Mesh velocities are EXACT at t(n+1)', 'Method 2: Mesh velocities are EXACT at t(n+1/2)','Method 3: Mesh velocities satisfy DGCL')
xlabel('dt [s]')
ylabel('Error 1: (RMS at final time)')

% error 2 plots 
figure(2);
hold("off");
loglog(dt_arr, e2_arr(:,1), 'b-x'); % method 1
hold("on");
loglog(dt_arr, e2_arr(:,2), 'r-x'); % method 2
loglog(dt_arr, e2_arr(:,3), 'g-x'); % method 3
legend('Method 1: Mesh velocities are EXACT at t(n+1)', 'Method 2: Mesh velocities are EXACT at t(n+1/2)','Method 3: Mesh velocities satisfy DGCL')
xlabel('dt [s]')
ylabel('Error 2: (RMS for entire simulation)')
