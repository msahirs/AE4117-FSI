clear all;

%% Setup

NdtP = 20; % number of steps per period
Np   = 0.25;  % number of periods to compute

plotsol = 1;    % (=1) plot mesh quality every time step
fmt     = '-b'; % linetype for plotting mesh quality in time

N = 11;  % Number of cells in y
M = 11;  % Number of cells in x

absdisp = 1; % (=1) absolute mesh deformation displacements

maxAngle = 45;  % Maximum rotation angle (degrees)


%% Standard settings
Lx = 1; % Length of domain in x
Ly = 1; % Length of domain in y

%% Initialization

% Compute node location
nodes = zeros((N+1)*(M+1),2);
for x=0:M
    nodes(x*(N+1)+1:(x+1)*(N+1),1) = x/M*Lx;
    nodes(x*(N+1)+1:(x+1)*(N+1),2) = [0:N]'/N*Ly;
end

% Compute mesh properties
cells        = create_mesh(N,M,nodes);        % Nodes per cell

% Mesh dynamics
omega = 2*pi;        % radial angular frequency
dt    = 1/NdtP;      % time step
rc    = [Lx/2 Ly/2]; % rotation center

nodes0 = nodes; % store mesh at t=0

% Moving nodes
id1     = (N+1)*floor((M-1)/2)+floor((N+1)/2);
id2     = id1+1;
id3     = id1+(N+1);
id4     = id3+1;
dispIDs = [id1 id2 id3 id4];

% Static node IDs
staticIDs = [1:N+1 (1:M-1)*(N+1)+1 (2:M)*(N+1) M*(N+1)+1:(M+1)*(N+1)];

% minimum mesh quality
qualOfMesh = ones(NdtP*Np+1,1);

%% Simulation
for iTimeStep=1:NdtP*Np
    
    % time and angular displacement
    t         = iTimeStep*dt;
    theta_old = pi/180*maxAngle*sin(omega*(t-dt));
    theta     = pi/180*maxAngle*sin(omega*t);
        
    % Mesh update
    displacements = zeros((N+1)*(M+1),3);
    for nodeID=dispIDs
        if (absdisp)
            xold   = nodes0(nodeID,1);
            yold   = nodes0(nodeID,2);
            dTheta = theta;
        else
            xold   = nodes(nodeID,1);
            yold   = nodes(nodeID,2);
            dTheta = theta-theta_old;
        end
        xnew = rc(1) + (xold-rc(1))*cos(dTheta) - (yold-rc(2))*sin(dTheta);
        ynew = rc(2) + (yold-rc(2))*cos(dTheta) + (xold-rc(1))*sin(dTheta);
        displacements(nodeID,:) = [1 xnew-xold ynew-yold];
    end
    for nodeID=staticIDs
        displacements(nodeID,:) = [1 0 0];
    end
    if (absdisp)
        nodes   = movemesh(nodes0,displacements);
    else
        nodes   = movemesh(nodes,displacements);
    end
    
    % determine error compared to uniform flow
    qualOfMesh(iTimeStep+1) = min(plotMesh(cells,nodes));

    if (plotsol)
        close(1);
        figure(1);
        plotMesh(cells,nodes);
        figure(1);
%        pause
    end

    
    
end

%% plot results

figure(2);
hold on;
title('Mesh Quality');
ylabel('min orthogonality');
xlabel('time');
plot((0:NdtP*Np)*dt,qualOfMesh,fmt);

figure(1);
close(1);
figure(1);
title('Mesh quality (orthogonality)');
xlabel('x');
ylabel('y');
colorbar;
plotMesh(cells,nodes);
