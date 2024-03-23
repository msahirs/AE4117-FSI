function FaceSn = compute_faceSn(cell, nodes)
    FaceSn = zeros(4,2);
    
    % Compute face surface area * face normal
    for i=1:4
        FaceSn(i,:) = ([[0 -1];[1 0]]*[nodes(cell(mod(i,4)+1),:)-nodes(cell(i),:)]')';
    end
end