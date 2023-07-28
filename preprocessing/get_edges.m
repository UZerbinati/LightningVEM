%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: get_edges 
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% Returns the indices of the edges that join an array of vertices 
%
% Input
% =====
% edges  : Edge matrix
% vertex : Vertices of the edges 
%
% Output
% ======
% index : Indices of the edge
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% May  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function [index] = get_edges(edges,vertex)

    nedges = numel(vertex);                                                                          %Number of edges
    index  = zeros(nedges,1);                                                                        %Initializate the index vector
    
    for i=1:(nedges-1)
        
        a = vertex(i);                                                                               %Index ot the i-th     vertex
        b = vertex(i+1);                                                                             %Index of the i-th + 1 vertex
        
        if (a > b)                                                                                   %If a is bigger than b switch
            
            temp = a;
            a    = b;
            b    = temp;
            
        end
        
        index(i,1) = edges(a,b);                                                                     %Get the name of the edge joining a and b
        
    end
    
        a = vertex(end);                                                                             %Same procedure for the first and the last vertex
        b = vertex(1);
        
        if (a > b)
            
            temp = a;
            a    = b;
            b    = temp;
            
        end
        
        index(end,1) = edges(a,b);
    
    return
end