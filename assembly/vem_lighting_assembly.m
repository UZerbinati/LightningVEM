function [K_global, M_global, f_global,domainMesh] = vem_lighting_assembly(domainMesh, matProps, f, k)

% vem_lighting_assembly: This function constructs the global matrix
% associated to the discretization of the Advection-Diffusion-Equation with
% the lighting virtual element method.
%
% Input parameters:
% domainMesh: the structure that contains the information on the mesh;
%   matProps: the list of parameters for the PDE;
%          f: the right-hand side of the PDE;
%          k: the degree of the polynomials (actaully = 1).
%
% Output parameters:
%   K_global: the global matrix;
%   f_global: the right-hand side of the linear system;
% domainMesh: the updated structure.

%% MEMORY ALLOCATION
siz = domainMesh.nvertex + (k - 1)*domainMesh.nedges + (k-1)*k/2 * domainMesh.npolygon;              %Size of the linear system     

K_global = zeros(siz,siz);                                                                           %Memory allocation
M_global = zeros(siz,siz);
f_global = zeros(siz, 1);                                                                          

domainMesh.polygon  = {1:domainMesh.npolygon};                                                       %Allocate memory for polygon information
domainMesh.diameter = 0;                                                                             %Initialize diameter of the mesh

connect  = domainMesh.connect;
coords   = parallel.pool.Constant(domainMesh.coords);
polygons = cell(domainMesh.npolygon,1);
K_loc    = cell(domainMesh.npolygon,1);
f_loc    = cell(domainMesh.npolygon,1);

%% CONSTRUCT THE LINEAR SYTEM
parfor i = 1:domainMesh.npolygon
   
    local_vertex = connect{i,1};                                                                     %Get the vertex indexes of the i-th polygon

    vertex = coords.Value(local_vertex,:);

    [K_loc{i}, M_loc{i}, f_loc{i}, polygons{i}] = vem_lighting_element(vertex, matProps, f, k);      %Assembly local stiffness matrix and load term

end

for i = 1:domainMesh.npolygon
    
    polygons{i}.local_dofs = get_local_dofs(domainMesh, k, i);                                       %Store the polygon information 
    local_dofs = polygons{i}.local_dofs;

    if (polygons{i}.diameter > domainMesh.diameter)                                                  %Update the diameter of the mesh

        domainMesh.diameter = polygons{i}.diameter;

    end

    K_global(local_dofs, local_dofs) = K_global(local_dofs, local_dofs) + K_loc{i};
    M_global(local_dofs, local_dofs) = M_global(local_dofs, local_dofs) + M_loc{i};
    f_global(local_dofs)             = f_global(local_dofs)             + f_loc{i};

end

domainMesh.polygon = polygons;