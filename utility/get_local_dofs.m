function [local_dofs] = get_local_dofs(domainMesh, polynomial, index, local_vertex, local_edges)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: get_local_dofs
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function creates an array that contains the degrees of freedom associated to the index-th 
% polygon
%
% Input
% =====
% domainMesh   : The struct relative to the mesh that we are using
% polynomial   : The struct that contains the information on the polynomials
% index        : The index of the polygon
% local_vertex : The array that contains the indexes of the vertices
% local_edges  : The array that contains the indexes of the edges
%
% Output
% ======
% local_dofs   : The array that contains the indexes of the local DOFs
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mai  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

k = polynomial.k;                                                                                    %Degree of the polynomials

nvertex = numel(local_vertex);                                                                       %Number of vertex
size    = nvertex * k + polynomial.int;                                                              %Size of the DOFs vector

local_dofs            = zeros(size,1);                                                               %Initialization            
local_dofs(1:nvertex) = local_vertex;                                                                %The first DOFs are the indexes of the vertices
local_vertex          = [local_vertex; local_vertex(1)];

if (k > 1)
    
    for i = 1:nvertex

        if (local_vertex(i) < local_vertex(i+1))
            local_dofs(nvertex+1+(i-1)*(k-1):nvertex+i*(k-1)) = (domainMesh.nvertex ...                   %The second DOFs are associated to the pointwise values
                                                              + 1 + (local_edges(i)-1)*(k-1):domainMesh.nvertex + local_edges(i)*(k-1));                  %on the interior of the edges
        else

            local_dofs(nvertex+1+(i-1)*(k-1):nvertex+i*(k-1)) = (domainMesh.nvertex + local_edges(i)*(k-1):-1: ...                   %The second DOFs are associated to the pointwise values
                                                              domainMesh.nvertex+ 1 + (local_edges(i)-1)*(k-1));                  %on the interior of the edges
        end


    end
    
    local_dofs(nvertex*k+1 : end) = domainMesh.nvertex + domainMesh.nedges * (k-1) ...               %The third DOFs are associated to the moments in the 
                                  + (index-1) * polynomial.int + (1:polynomial.int);                                          %interior of the domain
end

end

