function [Gs] = vem_matrix_lighting_Gs(grad_x, grad_y, p_quad) 

siz = size(grad_x,4);
Gs   = zeros(2*siz);

for i = 1:siz

    for j = 1:siz

        integrand = grad_x(:,:,:,i) .* grad_x(:,:,:,j) + 0.5 * grad_y(:,:,:,i) .* grad_y(:,:,:,j);
        Gs(i,j) = quadrature_2D(p_quad, integrand,"Evaluated");

        integrand = 0.5 * grad_x(:,:,:,i) .* grad_y(:,:,:,j);
        Gs(i,j+siz) = quadrature_2D(p_quad, integrand,"Evaluated");

        integrand = 0.5 * grad_y(:,:,:,i) .* grad_x(:,:,:,j);
        Gs(i+siz,j) = quadrature_2D(p_quad, integrand,"Evaluated");

        integrand = grad_y(:,:,:,i) .* grad_y(:,:,:,j) + 0.5 * grad_x(:,:,:,i) .* grad_x(:,:,:,j);
        Gs(i+siz,j+siz) = quadrature_2D(p_quad, integrand,"Evaluated");

    end

end

end