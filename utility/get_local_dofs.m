function [local_dofs, vertex, boundary_edges] = get_local_dofs(Mesh, k, index)

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

i_vertex = Mesh.connect{index,1};
i_edges  = get_edges(Mesh.edges, i_vertex);
dimm     = (k-1)*k/2;                                                                            %Degree of the polynomials

nv   = numel(i_vertex);                                                                              %Number of vertex
size = nv * k + dimm;                                                               %Size of the DOFs vector

local_dofs       = zeros(size,1);                                                                    %Initialization            
local_dofs(1:nv) = i_vertex;                                                                         %The first DOFs are the indexes of the vertices
i_vertex         = [i_vertex; i_vertex(1)];

boundary_edges = zeros(1, numel(i_edges));

if (k > 1)
    
    for i = 1:nv

        range = nv+1+(i-1)*(k-1):nv+i*(k-1);

        if (i_vertex(i) < i_vertex(i+1))
            
            local_dofs(range) = Mesh.nvertex + ((i_edges(i)-1)*(k-1)+1 : i_edges(i)*(k-1));                           
        
        else
        
            local_dofs(range) = Mesh.nvertex + ((i_edges(i)*(k-1) : -1 : (i_edges(i)-1)*(k-1)+1)); 
       
        end


    end
    
    local_dofs(nv*k+1 : end) = Mesh.nvertex + Mesh.nedges * (k-1) ...                                %The third DOFs are associated to the moments in the 
                             + (index-1) * dimm + (1:dimm);        %interior of the domain

end

for i = 1:numel(i_edges)

    boundary_edges(i) = ismember(i_edges(i), Mesh.boundary_edges);

end

vertex = Mesh.coords(i_vertex,:);

end