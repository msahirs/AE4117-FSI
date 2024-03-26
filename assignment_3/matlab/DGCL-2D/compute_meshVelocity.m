function meshVelocity = compute_meshVelocity(mesh,nodes,nodesTN,dt)
    N = length(mesh(:,1));  % number of cells in mesh
    
    % NOTE: Vectors are stored in row!
    
    if (0) % violate DGCL
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
            Sn_TN   = compute_faceSn(cell,nodesTN);   % face surfaces*normal at t_n
            c_TNP1  = compute_facecenter(cell,nodes);   % face centers at t_n+1
            Sn_TNP1 = compute_faceSn(cell,nodes);      % face surfaces*normal at t_n+1
            % c_TNP1 - c_TN

            for IDface=1:4
                %% IMPLEMENT THE DGCL HERE %%
                
                h_diff = c_TNP1(IDface,:) - c_TN(IDface,:);
                surf_avg = (Sn_TNP1(IDface,:) + Sn_TN(IDface,:))/2;
                % sweptVolume_dt = (h_diff.*surf_avg)/dt;  % element wise
                sweptVolume_dt = dot(h_diff,surf_avg)/dt;  % dot product
                velocity = sweptVolume_dt ./ Sn_TNP1(IDface,:) ;
                % velocity = (sweptVolume) ./ ((Sn_TNP1(IDface,:) - Sn_TN(IDface,:))*dt);
                
                if  (isnan(velocity(1)) || isinf(velocity(1)))
                    velocity(1) = 0;
                end
                if (isnan(velocity(2)) || isinf(velocity(2)))
                    velocity(2) = 0;
                end
            
                meshVelocity(IDcell,IDface,1) = velocity(1);
                meshVelocity(IDcell,IDface,2) = velocity(2);
            end
                       
        end
    end

end
