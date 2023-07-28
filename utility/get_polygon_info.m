function [polygon] = get_polygon_info(vertex) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: get_polygon_info
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function creates the struct polygon 
%
% Input
% =====
% vertex : The coordinates of the vertices
%
% Output
% ======
% polygon : A struct that contains the information on the polygon with the following fields:
%
%       area     : The area of the polyogon
%       diameter : The diameter of the polygon
%       centroid : The centroid of the polygon
%       edges    : The lengths of the edges
%       nedges   : The space dimension
%       paramx   : The parametrizations of the edges (x component)
%       paramy   : The parametrizations of the edges (y component)
%       vnormal  : The normal vectors of the polygon
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mai  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                                                                            
    polygon.nedges    = (numel(vertex) - 2) / 2;                                                     %Number of edges
    polygon.area      = 0;                                                                           %Area of the polygon
    polygon.edges     = zeros(polygon.nedges, 1);                                                    %Edges length
    polygon.vnormal   = zeros(polygon.nedges, 2);                                                    %Normal vectors
    polygon.centroid  = [0 0];                                                                       %Centroid
    polygon.diameter  = 0;                                                                           %Diameter of the polygon
   
    for i=1:polygon.nedges
        
        polygon.area          = polygon.area + vertex(i,1) * vertex(i+1,2) ...
                                             - vertex(i,2) * vertex(i+1,1);
                          
        polygon.edges(i)      = sqrt(sum( (vertex(i,:) - vertex(i+1,:)).^2 ));
         
        polygon.vnormal(i,:)  = [vertex(i+1,2) - vertex(i,2), -vertex(i+1,1) + vertex(i,1)];
                            
        polygon.vnormal(i,:)  = polygon.vnormal(i,:) ...
                              / hypot(polygon.vnormal(i,1),polygon.vnormal(i,2));
                          
        polygon.centroid(1)   = polygon.centroid(1) + (vertex(i,1)  + vertex(i+1,1)) ... 
                              * (vertex(i,1)  * vertex(i+1,2) - vertex(i+1,1) * vertex(i,2));
                          
        polygon.centroid(2)   = polygon.centroid(2) + (vertex(i,2)  + vertex(i+1,2)) ... 
                              * (vertex(i,1)  * vertex(i+1,2) - vertex(i+1,1) * vertex(i,2));
               
        for j = i+1:polygon.nedges
            
            d = sqrt(sum( (vertex(i,:) - vertex(j,:)).^2 ));
             
            if d > polygon.diameter
                
                polygon.diameter = d;
                 
            end
            
        end

    end
    
    polygon.perimeter = sum(polygon.edges);
    polygon.area = polygon.area/2;
    
    if (polygon.area<0)
        
        polygon.area = polygon.area * (-1);
        
    end
    
    polygon.centroid  = polygon.centroid / (6 * polygon.area);

end