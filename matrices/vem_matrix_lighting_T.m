function [T] = vem_matrix_lighting_T(val_int, grad_x, grad_y, b_val1, b_val2, p_quad)

dim = size(val_int,4);
T   = zeros(dim);

for i=1:dim

    for j=1:dim

        integrand = (b_val1 .* grad_x(:,:,:,j) + b_val2 .* grad_y(:,:,:,j)) .* val_int(:,:,:,i);

        T(i,j)    = quadrature_2D(p_quad, integrand,"Evaluated");

    end

end