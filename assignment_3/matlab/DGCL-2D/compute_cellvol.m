function vol = compute_cellvol(cells,nodes)
    N   = length(cells(:,1));
    vol = zeros(N,1);
    
    % Cell volume is summation of inner products of the face location
    % and its face surface * normal
    for iCell=1:N
        n1 = nodes(cells(iCell,1),:); % corner 1
        n2 = nodes(cells(iCell,2),:); % corner 2
        n3 = nodes(cells(iCell,3),:); % corner 3
        n4 = nodes(cells(iCell,4),:); % corner 4
        
        c = 0.25 * (n1 + n2 + n3 + n4); % cell center
        
        vol(iCell) = 0;
        
        FaceC  = compute_facecenter(cells(iCell,:), nodes); % face center
        FaceSn = compute_faceSn(cells(iCell,:), nodes);     % face surface * normal
        
        for iFace=1:4
            vol(iCell) = vol(iCell) + (FaceC(iFace,:) - c) * FaceSn(iFace,:)';
        end
        
        vol(iCell) = 0.5 * vol(iCell);
    end
end