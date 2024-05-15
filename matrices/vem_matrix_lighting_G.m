function [G] = vem_matrix_lighting_G(grad_x, grad_y, p_quad)

% vem_matrix_lighting_G: this computes the matrix that stores the integral grad(u_i)*grad(u_j).
%
% Input parameters:
% grad_x: x-derivative of the basis functions;
% grad_y: y-derivative of the basis functions;
% p_quad: quadrature points.
%
% Output parameters:
% G: the matrix.

siz = size(grad_x,4);
G   = zeros(siz);

for i = 1 : siz

    integrand = grad_x(:,:,:,i).^2 + grad_y(:,:,:,i).^2;
    G(i,i) = quadrature_2D(p_quad, integrand,"Evaluated");

    for j = i+1 : siz

        integrand = grad_x(:,:,:,i).*grad_x(:,:,:,j) + grad_y(:,:,:,i).*grad_y(:,:,:,j);
        G(i,j) = quadrature_2D(p_quad, integrand,"Evaluated");
        G(j,i) = G(i,j);

    end

end