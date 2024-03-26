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
                sweptVolume_dt = dot(h_diff,surf_avg)/dt;  % dot product
                % velocity = sweptVolume_dt ./ Sn_TNP1(IDface,:) ;
                
                % if  (isnan(velocity(1)) || isinf(velocity(1)))
                %     velocity(1) = 0;
                % end
                % if (isnan(velocity(2)) || isinf(velocity(2)))
                %     velocity(2) = 0;
                % end
            
                % meshVelocity(IDcell,IDface,1) = velocity(1);
                % meshVelocity(IDcell,IDface,2) = velocity(2);
                
                %% IMPLEMENT THE BONUS HERE %%
                c21 = 1767732205903 / 4055673282236;
                c22 = 1767732205903 / 4055673282236;
                c31 = 2746238789719 / 10658868560708;
                c32 = -640167445237 / 6845629431997;
                c33 = 1767732205903 / 4055673282236;
                c41  = 1471266399579 / 7840856788654; 
                c42 = -4482444167858 / 7529755066697;
                c43 = 11266239266428 / 11593286722821;
                c44 = 1767732205903 / 4055673282236;
                
                % Compute intermediate velocities
                velocity_intermediate = zeros(4, 1);
                velocity_intermediate(1) = sweptVolume_dt / Sn_TNP1(IDface,:);
                % velocity_intermediate(2) = (sweptVolume_dt - c21 * velocity_intermediate(1) - c22 * velocity_intermediate(2)) / Sn_TNP1(IDface,:);
                velocity_intermediate(2) = (sweptVolume_dt - c21 * velocity_intermediate(1)) / (Sn_TNP1(IDface,:) * (1+c22));
                velocity_intermediate(3) = (sweptVolume_dt - c31 * velocity_intermediate(1) - c32 * velocity_intermediate(2)) / (Sn_TNP1(IDface,:)*(1+c33));
                velocity_intermediate(4) = (sweptVolume_dt - c41 * velocity_intermediate(1) - c42 * velocity_intermediate(2) - c43 * velocity_intermediate(3)) / (Sn_TNP1(IDface,:)*(1+c44));
                
                % Compute final velocity
                velocity = velocity_intermediate(4);

                meshVelocity(IDcell,IDface,1) = velocity(1);
                meshVelocity(IDcell,IDface,2) = velocity(2);
            end                       
        end
    end    

end
