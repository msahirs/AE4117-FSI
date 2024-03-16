function cells = create_mesh(N,M,nodes)
    cells = zeros(N*M,4);
    
    iCell = 1;
    
    for j=1:M
        for i=1:N
            n1 = (j-1)*(N+1)+i;
            n2 = n1+1;
            n4 = j*(N+1)+i;
            n3 = n4 + 1;
            cells(iCell,:) = [n1 n2 n3 n4];
            iCell = iCell + 1;
        end
    end
end