function [boundary_vertices, boundary_intdofs] = get_boundary_dofs(domainMesh, boundary_vertex, k)

% get_boundary_dofs: obtain the indexes of the boundary dofs
%
% Input parameters:
%      domainMesh: the structure that contains the information on the mesh;
% boundary_vertex: indexes of the boundary vertices;
%               k: degree of the polynomials.
%
% Output parameters:
% boundary_vertices: the indexes of the first type of dofs;
%  boundary_intdofs: the indexes of the second type of dofs.

boundary_intdofs = zeros(numel(domainMesh.boundary_edges)*(k-1),1);                                  %Memory allocation for boundary DOFs

if (k > 1)

    for j = 1:numel(domainMesh.boundary_edges)
    
        i     = domainMesh.boundary_edges(j);

        boundary_intdofs(1 +(j-1)*(k-1) : j*(k-1)) = (domainMesh.nvertex + 1 + (i-1)*(k-1) : ...
                                                      domainMesh.nvertex     +  i   *(k-1));

    end

end

boundary_vertices =  [boundary_vertex; boundary_intdofs];

end