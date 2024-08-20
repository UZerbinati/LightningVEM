function [A1, A2, B1, B2, f_global,domainMesh] = vem_adr_assembly(domainMesh, matProps, polynomial, f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_adr_assembly
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function construct the linear system associated to the discretization of the 
% advection-diffusion-reaction problem and the associated right-hand side.
%
% Input
% =====
% domainMesh : The struct that contains the information on the mesh that we are using
% matProps   : The struct that contains the parameters of the equation
% polynomial : The struct that contains the information on the polynomials that we are using
% f          : The right-hand side of the linear system
%
% Output
% ======
% K_global   : The matrix of the linear system
% f_global   : The right-hand side
% domainMesh : The updated struct with polygon info
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MEMORY ALLOCATION
siz = domainMesh.nvertex + (polynomial.k - 1) * domainMesh.nedges ...                                %Size of the linear system
    + polynomial.int * domainMesh.npolygon;         

A1 = zeros(siz,siz);                                                                           %Memory allocation
A2 = zeros(siz,siz);
B1 = zeros(siz,siz);
B2 = zeros(siz,siz);

f_global = zeros(siz, 1);                                                                          

domainMesh.polygon  = {1:domainMesh.npolygon};                                                       %Allocate memory for polygon information
domainMesh.diameter = 0;                                                                             %Initialize diameter of the mesh

%% CONSTRUCT THE LINEAR SYTEM
for i=1:domainMesh.npolygon
    
    local_dofs   = get_local_dofs(domainMesh, polynomial.k, i);                                      %Get the DOFs indexes of the i-th polygon
    
    [A1_loc, A2_loc, B1_loc, B2_loc, f_loc, polygon] = vem_adr_element(domainMesh, matProps,...
                                                                       polynomial, i, f);            %Assembly local stiffness matrix and load term
                                                                       
    polygon.local_dofs    = local_dofs;                                                              %Store local_dofs in the struct polygon
    domainMesh.polygon{i} = polygon;                                                                 %Store the polygon information 
    
    if (domainMesh.polygon{i}.diameter > domainMesh.diameter)                                        %Update the diameter of the mesh

        domainMesh.diameter = domainMesh.polygon{i}.diameter;

    end
    
    A1(local_dofs, local_dofs) = A1(local_dofs, local_dofs) + A1_loc;                                %Update the matrix
    A2(local_dofs, local_dofs) = A2(local_dofs, local_dofs) + A2_loc;
    B1(local_dofs, local_dofs) = B1(local_dofs, local_dofs) + B1_loc;
    B2(local_dofs, local_dofs) = B2(local_dofs, local_dofs) + B2_loc;
    f_global(local_dofs)       = f_global(local_dofs)       + f_loc;                                 %Update the load term
end

end


