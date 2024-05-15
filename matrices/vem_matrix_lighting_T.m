function [T] = vem_matrix_lighting_T(val_int, grad_x, grad_y, b_val1, b_val2, p_quad)

% vem_matrix_lighting_T: this computes the matrix that stores the integral beta*grad(u_i)*u_j.
%
% Input parameters:
% val_int: pointwise values of the basis functions;
%  b_val1: poinwise values of the first component of beta;
%  b_val2: poinwise values of the second component of beta;
%  p_quad: quadrature points.
%
% Output parameters:
% T: the matrix.

dim = size(val_int,4);
T   = zeros(dim);

for i=1:dim

    for j=1:dim

        integrand = (b_val1 .* grad_x(:,:,:,j) + b_val2 .* grad_y(:,:,:,j)) .* val_int(:,:,:,i);

        T(i,j)    = quadrature_2D(p_quad, integrand,"Evaluated");

    end

end