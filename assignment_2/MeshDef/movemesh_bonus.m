function nodesNew = movemesh_bonus(nodes,disp,r)
    
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
        Mb = zeros(Nb, Nb);
        P = zeros(Nb,3);
        P = [ones(Nb, 1), Xb(:,1), Xb(:,2)];
        for i = 1:Nb
            for j = 1:Nb
                distance = sqrt((Xb(i,1) - Xb(j,1))^2 + (Xb(i,2) - Xb(j,2))^2);
                if distance > r
                    Mb(i,j) = 0;
                else
                    Mb(i,j) = (1-distance/r)^4*(4*distance/r+1);
                end
            end
        end
        
        rbfMat = [Mb P ; P' zeros(3)];

        Mb = zeros(Ni, Nb);
        P = zeros(Ni, 3);
        P = [ones(Ni, 1), Xi(:,1), Xi(:,2)];

        % evaluate RBF function for each internal node
        for i = 1:Ni
            for j = 1:Nb
                distance = sqrt((Xi(i,1) - Xb(j,1))^2 + (Xi(i,2) - Xb(j,2))^2);
                if distance > r
                    Mb(i,j) = 0;
                else
                    Mb(i,j) = (1-distance/r)^4*(4*distance/r+1);
                end
            end
        end
        
        H = [Mb P]/rbfMat;
        H = H(:, 1:Nb);
        Di = H * Db;
        % Convert internal displacements to global disp vector
        for i=1:Ni
            disp(ID(i),2:3) = Di(i,:);
        end
        
        % Update node locations
        nodesNew = nodes+disp(:,2:3);
    end
end