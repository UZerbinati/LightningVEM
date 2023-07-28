 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  VEMLab
%         Source code  : http://camlab.cl/research/software/vemlab/
%           (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:               sort_boundary_edges_left 
%
% Created by : M. Trezzi
% Updated by :
% Updated by :
%
%--------------------------------------------------------------------------
% Purpose
% =======
% Sorts the vertices on the left boundary. The mesh generator returns the 
% left nodes ordered by indices. We want order them by y-coordinate
%
% Usage
% =====
% left = sort_boundary_edges_left(domainMesh.coords, 
%                                   domainMesh.boundary_nodes.left);
%
% Input
% =====
% coords : coordinates of the vertices
% x      : indices of the vertices
%
% Output
% ======
% x : ordered array
%
%--------------------------------------------------------------------------
% References 
% ==========
%
%--------------------------------------------------------------------------
% Function's updates history
% ==========================
% Apr.  4, 2021: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [x] = sort_boundary_edges_left(coords, x)

    i = 1;
        
    while (i < numel(x))
            
        if (coords(x(i),2) > coords(x(i+1),2))
            
            temp = x(i);
            x(i) = x(i+1);
            x(i+1) = temp;
            
            if (i == 1)
                i = 0;
            else
                i = i - 2;
            end
            
        end

        i = i + 1;
        
    end
    
return
end