function [Bb] = vem_matrix_Bb(beta_val1, beta_val2, base_val_int, grad_val_int, p_quad, polynomial)

dim = polynomial.dim;                                                                                %Memory allocation
Bb  = zeros(dim, dim);

for i=1:dim

    for j=1:dim

        integrand =  (grad_val_int(:,:,:,j).*beta_val1 + grad_val_int(:,:,:,j+dim).*beta_val2) ...
                  .* base_val_int(:,:,:,i);
        Bb(i,j)   = quadrature_2D(p_quad, integrand, "Evaluated");

    end

end