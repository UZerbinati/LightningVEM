function [b_quad] = boundary_quadrature(polygon, n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: boundary_quadrature
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function constructs a quadrature on a boundary that is exact of order 2n-1 for a polygon
%
% Input
% =====
% polygon : The polygon
% n       : Number of nodes for the 1D Gauss quadrature
%
% Output
% ======
% b_quad : A struct that contains the following fields
%        
%        x : x-coordinates
%        y : y-coordinates
%        w : The weights of the quadrature
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


vertex   = polygon.vertex;
nedges   = polygon.nedges;

b_quad.x = zeros(nedges,n);                                                                          %Memory allocation
b_quad.y = zeros(nedges,n);
tau      = gauss_points_1d(n);
b_quad.w = gauss_weights_1d(n);

for i=1:nedges
    
    x_i = @(t) ( vertex(i+1,1) - vertex(i,1) ) / 2 .*t + ( vertex(i+1,1) + vertex(i,1) ) / 2;        %Parametrization of the edges
    y_i = @(t) ( vertex(i+1,2) - vertex(i,2) ) / 2 .*t + ( vertex(i+1,2) + vertex(i,2) ) / 2;
    
    b_quad.x(i,:) = x_i(tau);                                                                        %Points on the edges
    b_quad.y(i,:) = y_i(tau);
         
end

end