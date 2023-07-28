function [boundary_dofs, boundary_intdofs] = get_boundary_dofs(domainMesh, boundary_vertex, k)

boundary_intdofs = zeros(numel(domainMesh.boundary_edges)*(k-1),1);                                  %Memory allocation for boundary DOFs

if (k > 1)

    for j = 1:numel(domainMesh.boundary_edges)
    
        i     = domainMesh.boundary_edges(j);
        index = numel(boundary_vertex) + 1 +(j-1)*(k-1):numel(boundary_vertex) + j*(k-1);

        boundary_intdofs(1 +(j-1)*(k-1) : j*(k-1)) = (domainMesh.nvertex + 1 + (i-1)*(k-1) : ...
                                                      domainMesh.nvertex     +  i   *(k-1));

    end

end

boundary_dofs =  [boundary_vertex; boundary_intdofs];

end