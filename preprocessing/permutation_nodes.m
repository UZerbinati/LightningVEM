function domainMesh = permutation_nodes(domainMesh, perm)

    nvertex   = domainMesh.nvertex;
    vertex    = [1:nvertex];
    intvertex = setdiff([1:nvertex], domainMesh.boundary_nodes.all);
    
    for i = 1:numel(intvertex)
        
        rand_x = 2 * perm * (rand - 0.5);
        rand_y = 2 * perm * (rand - 0.5);
        
        domainMesh.coords(intvertex(i),:) = domainMesh.coords(intvertex(i),:) + [rand_x rand_y];
        
    end
    
end