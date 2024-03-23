function connectivity = compute_connectivity(mesh)
    N = length(mesh(:,1));
 
    connectivity = zeros(N,4);
 
    % For each cell face, determine neighboring cells
    for iCell = 1:N
        tempmesh          = mesh;
        tempmesh(iCell,:) = 0;
        for iFace = 1:4
            id1 = mesh(iCell,iFace);
            id2 = mesh(iCell,mod(iFace,4)+1);
            connectivity(iCell,iFace) = [1:N] * (sum(tempmesh == id1,2) & sum(tempmesh == id2,2));
        end
    end
end
                