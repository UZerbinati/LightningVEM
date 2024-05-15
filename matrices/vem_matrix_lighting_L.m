function [L] = vem_matrix_lighting_L(grad_x, grad_y, p_quad) 

siz = size(grad_x,4);
L   = zeros(2*siz);

for i = 1:siz

    for j = 1:siz

        integrand = grad_x(:,:,:,i) .* grad_x(:,:,:,j);
        L(i,j) = quadrature_2D(p_quad, integrand,"Evaluated");

        integrand = grad_x(:,:,:,j) .* grad_y(:,:,:,i);
        L(i,j+siz) = quadrature_2D(p_quad, integrand,"Evaluated");
        
        integrand = grad_y(:,:,:,j) .* grad_x(:,:,:,i);
        L(i+siz,j) = quadrature_2D(p_quad, integrand,"Evaluated");

        integrand = grad_y(:,:,:,i) .* grad_y(:,:,:,j);
        L(i+siz,j+siz) = quadrature_2D(p_quad, integrand,"Evaluated");
            

        
    end

end

end