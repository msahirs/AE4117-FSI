function FaceC = compute_facecenter(cell, nodes)
    FaceC = zeros(4,2);
    
    % Compute center point of face
    for i=1:4
        FaceC(i,:) = 0.5 * (nodes(cell(i),:)+nodes(cell(mod(i,4)+1),:));
    end
end
    