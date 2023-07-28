 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  VEMLab
%         Source code  : http://camlab.cl/research/software/vemlab/
%           (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:              sort_boundary_edges_bottom 
%
% Created by : M. Trezzi
% Updated by :
% Updated by :
%
%--------------------------------------------------------------------------
% Purpose
% =======
% Sorts the vertices on the bottom boundary. The mesh generator returns the 
% bottom nodes ordered by indices. We want order them by x-coordinate
%
% Usage
% =====
% bottom = sort_boundary_edges_bottom(domainMesh.coords, 
%                                     domainMesh.boundary_nodes.bottom);
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

function [x] = sort_boundary_edges_bottom(coords, x)

    i = 1;
        
    while (i < numel(x))
            
        if (coords(x(i),1) > coords(x(i+1),1))

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