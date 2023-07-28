function [P0] = projection_P0(base_val_bound, polynomial, polygon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: projection_P0
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function computes the integral mean on the boundary of the polynomial basis function
%
% Input
% ===== 
% base_val_bound : Values of the basis functions on the boundary
% polynomial     : Information on the polynomial used to compute G
% polygon        : Information on the polygon
%
% Output
% ======
% P0 : The integral mean
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% May.  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

P0      = zeros(1,polynomial.dim);                                                                   %Memory allocation
weights = polynomial.weights;

for k = 1:polynomial.dim

    P0(k) = P0(k) + weights * base_val_bound(2:end,:,k)' * polygon.edges;

end

P0 = P0 ./ polygon.perimeter;

end