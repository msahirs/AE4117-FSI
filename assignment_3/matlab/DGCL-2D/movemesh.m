function nodesNew = movemesh(nodes,disp)
    
    Nb = sum(disp(:,1));       % boundary points (constraints)
    Ni = length(disp(:,1))-Nb; % internal points
    
    if (Ni == 0)
        % All points have prescribed displacements
        % Update node locations
        nodesNew = nodes+disp(:,2:3);
    else
        nodesNew = nodes+disp(:,2:3);
    end
end
