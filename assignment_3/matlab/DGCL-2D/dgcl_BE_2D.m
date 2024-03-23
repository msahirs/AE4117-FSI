clear all;

%% Setup

NdtP = 10; % number of steps per period
Np   = 2;  % number of periods to compute

plotsol = 1;    % (=1) plot solution every time step
fmt     = '-b'; % linetype for plotting error in uniform solution


%% Standard settings
Lx = 1; % Length of domain in x
Ly = 1; % Length of domain in y

N = 3;  % Number of cells in y
M = 3;  % Number of cells in x

absdisp = 1; % (=1) absolute mesh deformation displacements

%% Initialization

% Compute node location
nodes = zeros((N+1)*(M+1),2);
for x=0:M
    nodes(x*(N+1)+1:(x+1)*(N+1),1) = x/M*Lx;
    nodes(x*(N+1)+1:(x+1)*(N+1),2) = [0:N]'/N*Ly;
end

% Compute mesh properties
cells        = create_mesh(N,M,nodes);        % Nodes per cell
connectivity = compute_connectivity(cells);   % Connectivity per face
volumes      = compute_cellvol(cells,nodes);  % Cell volumes
meshVelocity = zeros([size(cells) 2]);        % mesh velocity

% Initial solution (uniform)
solution = ones(N*M,1);

% Show initial solution
figure(1);
plotSolution(cells,nodes,solution);

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

% maximum error compared to uniform solution
errInSol = zeros(NdtP*Np+1,1);

%% Simulation
for iTimeStep=1:NdtP*Np
    
    % time and angular displacement
    t         = iTimeStep*dt;
    theta_old = pi/180*45*sin(omega*(t-dt));
    theta     = pi/180*45*sin(omega*t);
    
    % store data at tn
    volTN   = volumes;
    nodesTN = nodes;
    
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
    volumes = compute_cellvol(cells,nodes);
    
    % Mesh velocities
    meshVelocity = compute_meshVelocity(cells,nodes,nodesTN,dt);
    meshVel      = zeros(1,2);
        
    % Setting up system to solve; for each cell:
    %
    %   (U*vol)^(n+1) - (U*vol)^n             /                 \(n+1)
    %   ------------------------- - SUM_faces | U * meshVel . Sn |     = 0
    %             dt                          \                 /face
    %
    % which can be written for the total system as:
    %
    % L U^(n+1) = R U^n
    %
    % We use a simple avarage to compute the solution at a cell face:
    %   U_face = 0.5 * (U_leftcell + U_rightcell)
    % 
    % If there is no connectivity (boundary faces), we just use the 
    % solution of the cell:
    %   U_face = U_cell
    
    L = zeros(N*M);
    R = zeros(N*M);
    
    for cellID=1:N*M
        L(cellID,cellID)  = volumes(cellID) / dt;
        R(cellID,cellID)  = volTN(cellID) / dt;
        faceSurfaceNormal = compute_faceSn(cells(cellID,:),nodes);
        for faceID=1:4
            faceSn         = faceSurfaceNormal(faceID,:);
            meshVel(1)     = meshVelocity(cellID,faceID,1);
            meshVel(2)     = meshVelocity(cellID,faceID,2);
            neighborCellID = connectivity(cellID,faceID);
            if (neighborCellID == 0)
                L(cellID,cellID) = L(cellID,cellID) - meshVel * faceSn';
            else
                L(cellID,cellID)         = L(cellID,cellID)         - 0.5 * meshVel * faceSn';
                L(cellID,neighborCellID) = L(cellID,neighborCellID) - 0.5 * meshVel * faceSn';
            end
        end
    end
    
    % solve for new solution at t_n+1
    solution = L\(R*solution);
    
    % determine error compared to uniform flow
    errInSol(iTimeStep+1) = max(abs(solution-1));

    if (plotsol)
        close(1);
        figure(1);
        plotSolution(cells,nodes,solution);
        figure(1);
%        pause
    end

    
    
end

%% plot results
close(1);
figure(1);
title('Solution');
ylabel('y');
xlabel('x');
colorbar;
plotSolution(cells,nodes,solution);

figure(2);
hold on;
title('Error');
ylabel('max |error|');
xlabel('time');
plot((0:NdtP*Np)*dt,errInSol,fmt);

figure(3);
close(3);
figure(3);
title('Mesh quality (orthogonality)');
xlabel('x');
ylabel('y');
colorbar;
plotMesh(cells,nodes);
