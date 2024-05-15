function [M] = vem_matrix_lighting_M(vem_val_int, p_quad)

% vem_matrix_lighting_M: this computes the matrix that stores the integral u_i*u_j.
%
% Input parameters:
% vem_val_int: pointwise values of the basis functions;
%      p_quad: quadrature points
%
% Output parameters:
% M: the matrix.

dim = size(vem_val_int,4);
M   = zeros(dim);

for i = 1:dim

    integrand = vem_val_int(:,:,:,i).^2;
    M(i,i)    = quadrature_2D(p_quad, integrand,"Evaluated");

    for j =  i+1:dim

        integrand = vem_val_int(:,:,:,i).*vem_val_int(:,:,:,j);
        M(i,j)    = quadrature_2D(p_quad, integrand,"Evaluated");
        M(j,i)    = M(i,j);

    end
    
end

end