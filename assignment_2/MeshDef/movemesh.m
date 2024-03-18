function nodesNew = movemesh(nodes,disp)
    
    Nb = sum(disp(:,1));       % boundary points (constraints)
    Ni = length(disp(:,1))-Nb; % internal points
    
    if (Ni == 0)
        % All points have prescribed displacements
        % Update node locations
        nodesNew = nodes+disp(:,2:3);
    else
        Xb = zeros(Nb,2); % coordinates boundary points
        Db = zeros(Nb,2); % displacements boundary points
        Xi = zeros(Ni,2); % coordinates internal points
        Di = zeros(Ni,2); % displacements internal points
        ID = zeros(Ni,1); % global ID for internal points
        ib = 1;
        ii = 1;
        
        % Split nodes and disp into boundary and internal points
        for i=1:Nb+Ni
            if (disp(i,1) == 1)
                Xb(ib,:) = nodes(i,:);
                Db(ib,:) = disp(i,2:3);
                ib       = ib+1;
            else
                Xi(ii,:) = nodes(i,:);
                ID(ii)   = i;
                ii       = ii+1;
            end
        end
        
        % Set up RBF interpolation matrix for constraints
        rbfMat = zeros(Nb+3);

        %% IMPLEMENT YOUR RBF MESH DEFORMATION HERE
        
        % Set up the RBF matrix and constraints
        % rbfMat = [ Mb P
        %            PT 0 ]        
        % constraints = displacements of the boundary nodes (separate in x
        % and y)        
        % determine rbf coefficients 
        % evaluate RBF function for each internal node
        % Convert internal displacements to global disp vector

        for i = 1:Nb
            for j = 1:Nb
                distance = sqrt((Xb(i,1) - Xb(j,1))^2 + (Xb(i,2) - Xb(j,2))^2);
                if distance == 0
                    rbfMat(i,j) = 0;
                else
                    rbfMat(i,j) = distance^2 * log(distance);
                end
            end
        end

        P = [ones(Nb, 1), Xb(:,1), Xb(:,2)];
        rbfMat(Nb+1:Nb+3, 1:Nb) = P';
        rbfMat(1:Nb, Nb+1:Nb+3) = P;

        constraints = [Db(:,1), Db(:,2); zeros(3, 2)];
        rbfCoeff = rbfMat \ constraints;

        % evaluate RBF function for each internal node
        for i = 1:Ni
            for j = 1:Nb
                distance = sqrt((Xi(i,1) - Xb(j,1))^2 + (Xi(i,2) - Xb(j,2))^2);
                if distance == 0
                    Di(i,:) = Di(i,:) + rbfCoeff(j, :) * 0;
                else
                    Di(i,:) = Di(i,:) + rbfCoeff(j) * distance^2 * log(distance);
                end
            end
        end

        % Convert internal displacements to global disp vector
        for i=1:Ni
            disp(ID(i),2:3) = Di(i,:);
        end
        
        % Update node locations
        nodesNew = nodes+disp(:,2:3);
    end
end