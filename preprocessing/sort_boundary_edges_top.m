 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  VEMLab
%         Source code  : http://camlab.cl/research/software/vemlab/
%           (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                sort_boundary_edges_top 
%
% Created by : M. Trezzi
% Updated by :
% Updated by :
%
%--------------------------------------------------------------------------
% Purpose
% =======
% Sorts the vertices on the top boundary. The mesh generator returns the 
% top nodes ordered by indices. We want order them by x-coordinate
%
% Usage
% =====
% top = sort_boundary_edges_top(domainMesh.coords, 
%                               domainMesh.boundary_nodes.top);
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

function [x] = sort_boundary_edges_top(coords, x)
    
    i = 1;                                                                  %Array index
        
    while (i < numel(x))                                                        
            
        if (coords(x(i),1) > coords(x(i+1),1))

            temp = x(i);                                                    %Switch
            x(i) = x(i+1);
            x(i+1) = temp;
            
            if (i == 1)                                                     %We want to check if the order is preserved. We decrease
                i = 0;                                                      %the index to check it.
            else
                i = i - 2;
            end
            
        end

        i = i + 1;
        
    end

return
end