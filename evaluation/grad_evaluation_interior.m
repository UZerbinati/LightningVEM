function [grad_val] = grad_evalutation_interior(base_val, p_quad, polynomial, polygon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: grad_evalutation_interior
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function evaluates the gradient of the basis function on the quadrature nodes in the 
% interior of a polygon
%
% Input
% ===== 
% grad       : The gradient of  the basis function
% p_quad     : The quadrature formula
% polynomial : The struct that contains the information on the polynomials that we are using
%
% Output
% ======
% grad_val : The pointwise values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

grad_val                           = zeros(size(p_quad.xi));
grad_val(:,:,:,1:2*polynomial.dim) = 0;

diffx = p_quad.xi  - polygon.centroid(1);
diffy = p_quad.eta - polygon.centroid(2);

%grad_val(:,:,:,1:polynomial.dim)     = base_val ./ diffx .* reshape(polynomial.deg(:,1),1,1,1,[]);
%grad_val(:,:,:,polynomial.dim+1:end) = base_val ./ diffy .* reshape(polynomial.deg(:,2),1,1,1,[]);

for j = 2:polynomial.dim
    [derx, dery] = grad_fun(p_quad.xi,p_quad.eta,polynomial.deg(j,:),polygon);
    grad_val(:,:,:,j+polynomial.dim) = dery;
    grad_val(:,:,:,j)                = derx;
end

end