function [boundary_edges] = get_boundary_edges(domainMesh)

% get_boundary_edges: returns the indices of the edges on the boundary 
%
% Input parameters
% domainMesh: mesh of which we want the boundary edges
%
% Output parameters
% boundary_edges: array containing the indices of the boundary edges

    top    = sort_boundary_edges_top(domainMesh.coords, domainMesh.boundary_nodes.top);              %Sort the vertices on the top boundary
                                     
    bottom = sort_boundary_edges_bottom(domainMesh.coords, domainMesh.boundary_nodes.bottom);        %Sort the vertices on the bottom boundary
                                        
    right  = sort_boundary_edges_right(domainMesh.coords, domainMesh.boundary_nodes.right);          %Sort the vertices on the right boundary
                                       
    left   = sort_boundary_edges_left(domainMesh.coords, domainMesh.boundary_nodes.left);            %Sort the vertices on the left boundary
                                      
    
    top_boundary    = get_edges(domainMesh.edges, top)';                                             %Get edges indices
    bottom_boundary = get_edges(domainMesh.edges, bottom)';
    left_boundary   = get_edges(domainMesh.edges, left)';
    right_boundary  = get_edges(domainMesh.edges, right)';
    
    boundary_edges = [top_boundary(1:end-1)  bottom_boundary(1:end-1) ...                            %Assembly the result array
                      left_boundary(1:end-1) right_boundary(1:end-1)]';
    
end