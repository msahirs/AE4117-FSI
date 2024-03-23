function plotMesh(mesh,nodes)
    hold off;
    
    N = length(mesh(:,1));
    
    minC = 0.95;
    maxC = 0;
    
    for iCell=1:N
        Sn = compute_faceSn(mesh(iCell,:),nodes);
        n  = nodes(mesh(iCell,:),:);
        q  = 1;
        for iFace=1:4
            v1    = Sn(iFace,:);
            v2    = n(mod(iFace+1,4)+1,:) - n(mod(iFace,4)+1,:);
            angle = (v1*v2') / sqrt(v1*v1'*v2*v2');
            q     = min(q,(1-angle)/2);
        end
        minC = min(minC,q);
        maxC = max(maxC,q);
        X = [ [n(1,1) n(2,1)] ; [n(4,1) n(3,1)] ];
        Y = [ [n(1,2) n(2,2)] ; [n(4,2) n(3,2)] ];
        C = ones(2)*q;
        surface(X,Y,C);
        hold on;
    end
    caxis([minC maxC]);
end