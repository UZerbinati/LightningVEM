function [grad_val] = grad_evalutation_boundary(base_val, b_quad, polynomial, polygon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: grad_evalutation_boundary
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function evaluates the gradient of the basis on the quadrature nodes on the boundary 
% of a polygon
%
% Input
% ===== 
% grad       : The basis function
% b_quad     : The quadrature formula
% polynomial : The struct that contains the information on the polynomials that we are using
%
% Output
% ======
% grad_val : The pointwise values
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Apr.  4, 2021: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nedges   = size(b_quad.x,1);                                                                         %Number of edges
grad_val = zeros(size(b_quad.x));                                                                    %Memory allocation
grad_val = [grad_val(1:2,:); grad_val; grad_val];                                                    %Duplicate and add two rows
grad_val(:,:,1:polynomial.dim) = 0;                                                                  %Add a dimension for each row

diffx = b_quad.x - polygon.centroid(1);
diffy = b_quad.y - polygon.centroid(2);

%grad_val(2:nedges+1,:,:)   = base_val(2:end,:,:) ./ diffx .* reshape(polynomial.deg(:,1),1,1,[]);
%grad_val(3+nedges:end,:,:) = base_val(2:end,:,:) ./ diffy .* reshape(polynomial.deg(:,2),1,1,[]);

for j = 2:polynomial.dim
    [derx, dery] = grad_fun(b_quad.x,b_quad.y,polynomial.deg(j,:),polygon);
    grad_val(2:nedges+1,:,j)   = derx;
    grad_val(3+nedges:end,:,j) = dery;
end

grad_val(1,:)          = grad_val(nedges+1,:) ;
grad_val(2+nedges,:)   = grad_val(end,:);

end