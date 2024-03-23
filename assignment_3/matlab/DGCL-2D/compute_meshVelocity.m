function meshVelocity = compute_meshVelocity(mesh,nodes,nodesTN,dt)
    N = length(mesh(:,1));  % number of cells in mesh
    
    % NOTE: Vectors are stored in row!
    
    if (1) % violate DGCL
        % Base mesh velocities on displacement of face centers
        for IDcell=1:N
            cell    = mesh(IDcell,:);                   % cell definition
            c_TN    = compute_facecenter(cell,nodesTN); % face centers at t_n
            c_TNP1  = compute_facecenter(cell,nodes);   % face centers at t_n+1
            for IDface=1:4
                velocity = (c_TNP1(IDface,:) - c_TN(IDface,:)) / dt;
                meshVelocity(IDcell,IDface,1) = velocity(1);
                meshVelocity(IDcell,IDface,2) = velocity(2);
            end
        end
    else % use DGCL
        % Base mesh velocities on swept volumes of mesh faces
        for IDcell=1:N
            cell    = mesh(IDcell,:);                   % cell definition
            c_TN    = compute_facecenter(cell,nodesTN); % face centers at t_n
            Sn_TN   = compute_faceSn(cell,nodesTN);     % face surfaces*normal at t_n
            c_TNP1  = compute_facecenter(cell,nodes);   % face centers at t_n+1
            Sn_TNP1 = compute_faceSn(cell,nodes);       % face surfaces*normal at t_n+1
            for IDface=1:4

                %% IMPLEMENT THE DGCL HERE %%
            
            end
        end
    end

end
