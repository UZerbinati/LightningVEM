function [base_val] = base_evalutation_interior(p_quad, polynomial, polygon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: base_evalutation_interior
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function evaluates the basis function on the quadrature nodes in the  interior of a polygon
%
% Input
% ===== 
% base       : The basis function
% p_quad     : The quadrature formula
% polynomial : The struct that contains the information on the polynomials
%
% Output
% ======
% base_val : The pointwise values
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% May  07, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

base_val = ones(size(p_quad.xi));

base_val(:,:,:,size(polynomial.deg,1)) = 1;

% for d = 1:polynomial.dim
% 
%     base_val(:,:,:,d) = base(p_quad.xi,p_quad.eta,polynomial.deg(d,:));
%     
% end

diffx = (p_quad.xi  - polygon.centroid(1)) ./ polygon.diameter;
diffy = (p_quad.eta - polygon.centroid(2)) ./ polygon.diameter;

for d = 2:polynomial.dim
    if polynomial.deg(d,1) ~= 0
        base_val(:,:,:,d) = base_val(:,:,:,d-sum(polynomial.deg(d,:))) .* diffx;
    else
        base_val(:,:,:,d) = base_val(:,:,:,d-sum(polynomial.deg(d,:))-1) .* diffy;
    end
end

end