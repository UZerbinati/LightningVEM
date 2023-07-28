function [K_global,f_global,domainMesh] = vem_diffusion_reaction_assembly(domainMesh, matProps, polynomial, f, g, opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_diffusion_convection_assembly
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function construct the linear system associated to the discretization of the problem and 
% also the load term
%
% Input
% =====
% domainMesh : The struct that contains the information on the mesh that we are using
% matProps   : The struct that contains the parameters of the equation
% polynomial : The struct that contains the information on the polynomials that we are using
% f          : The right-hand side of the linear system
% opt        : Variable that contains the options of the method
%
% Output
% ======
% K_global   : The matrix of the linear system
% f_global   : The right-hand side
% domainMesh : The updated struct with polygon info
%
% Function's updates history
% ==========================
% Jul  20, 2022: insert opt
% May  07, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MEMORY ALLOCATION
siz = domainMesh.nvertex + (polynomial.k - 1) * domainMesh.nedges ...                                %Size of the linear system
    + polynomial.int * domainMesh.npolygon;         

K_global = zeros(siz,siz);                                                                           %Memory allocation
f_global = zeros(siz, 1);                                                                          

domainMesh.polygon  = {1:domainMesh.npolygon};                                                       %Allocate memory for polygon information
domainMesh.diameter = 0;                                                                             %Initialize diameter of the mesh

%% CONSTRUCT THE LINEAR SYTEM
for i=1:domainMesh.npolygon
    
    
    local_vertex = domainMesh.connect{i,1};                                                          %Get the vertex indexes of the i-th polygon
    local_edges  = get_edges(domainMesh.edges, domainMesh.connect{i,1});                             %Get the edge indexes of the i-th polygon
    local_dofs   = get_local_dofs(domainMesh, polynomial, i, local_vertex, local_edges);             %Get the DOFs indexes of the i-th polygon
    
    [K_loc, f_loc, polygon] = vem_diffusion_reaction_element(domainMesh, matProps, ...               %Assembly local stiffness matrix and load term
                                                             polynomial, i, f);          
         
    polygon.local_dofs    = local_dofs;                                                              %Store local_dofs in the struct polygon
    domainMesh.polygon{i} = polygon;                                                                 %Store the polygon information 
    
    if (domainMesh.polygon{i}.diameter > domainMesh.diameter)                                        %Update the diameter of the mesh

        domainMesh.diameter = domainMesh.polygon{i}.diameter;

    end
    
    K_global(local_dofs, local_dofs) = K_global(local_dofs, local_dofs) + K_loc;                     %Update the matrix
    f_global(local_dofs)             = f_global(local_dofs)             + f_loc;                     %Update the load term
end

end


