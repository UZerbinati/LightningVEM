function [K_global, M_global, f_global,domainMesh] = elasticity_assembly_cr(domainMesh, matProps, f1, f2, k)

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
siz = 2 * domainMesh.nvedges ;                                                                  %Size of the linear system     

K_global = zeros(siz,siz);                                                                           %Memory allocation
M_global = zeros(siz,siz);
f_global = zeros(siz, 1);                                                                          

domainMesh.polygon  = {1:domainMesh.npolygon};                                                       %Allocate memory for polygon information
domainMesh.diameter = 0;                                                                             %Initialize diameter of the mesh

connect  = domainMesh.connect;
coords   = parallel.pool.Constant(domainMesh.coords);
polygons = cell(domainMesh.npolygon,1);
K_loc    = cell(domainMesh.npolygon,1);
M_loc    = cell(domainMesh.npolygon,1);
f_loc    = cell(domainMesh.npolygon,1);

%% CONSTRUCT THE LINEAR SYTEM
for i = 1:domainMesh.npolygon
   
    local_vertex = connect{i,1};                                                                     %Get the vertex indexes of the i-th polygon

    vertex = coords.Value(local_vertex,:);

    [K_loc{i}, M_loc{i}, f_loc{i}, polygons{i}] = elasticity_element_cr(vertex, matProps, f1 ,f2, k);                %Assembly local stiffness matrix and load term

end

for i = 1:domainMesh.npolygon
    
    local_vertex = connect{i,1};                                                                     %Indexes of the vertices of the i-th polygon
    i_edges      = get_edges(Mesh.edges, local_vertex);

    local_dofs = [i_edges; i_edges + domainMesh.nedges];
    polygons{i}.local_dofs = local_dofs;                                                           %Store the polygon information 
    
    if (polygons{i}.diameter > domainMesh.diameter)                                                  %Update the diameter of the mesh

        domainMesh.diameter = polygons{i}.diameter;

    end

    K_global(local_dofs, local_dofs) = K_global(local_dofs, local_dofs) + K_loc{i};
    M_global(local_dofs, local_dofs) = M_global(local_dofs, local_dofs) + M_loc{i};
    f_global(local_dofs)             = f_global(local_dofs)             + f_loc{i};

end

domainMesh.polygon = polygons;