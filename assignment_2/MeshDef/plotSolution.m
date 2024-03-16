function plotSolution(mesh,nodes,solution)
    hold off;
    
    N = length(mesh(:,1));
    
    CAXIS([0.8 1.2]);
    for iCell=1:N
        n = nodes(mesh(iCell,:),:);
        X = [ [n(1,1) n(2,1)] ; [n(4,1) n(3,1)] ];
        Y = [ [n(1,2) n(2,2)] ; [n(4,2) n(3,2)] ];
        C = ones(2)*solution(iCell);
        surface(X,Y,C);
        hold on;
    end
end