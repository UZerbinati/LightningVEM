function [domainMesh] = fix_boundary_points(domainMesh)

    for i=1:numel(domainMesh.boundary_nodes.top)
        domainMesh.coords(domainMesh.boundary_nodes.top(i),2) = 1;
    end
    
    for i=1:numel(domainMesh.boundary_nodes.bottom)
        domainMesh.coords(domainMesh.boundary_nodes.bottom(i),2) = 0;
    end
    
    for i=1:numel(domainMesh.boundary_nodes.right)
        domainMesh.coords(domainMesh.boundary_nodes.right(i),1) = 1;
    end
    
    for i=1:numel(domainMesh.boundary_nodes.left)
        domainMesh.coords(domainMesh.boundary_nodes.left(i),1) = 0;
    end

end