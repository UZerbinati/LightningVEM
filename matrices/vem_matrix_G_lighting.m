function [G] = vem_matrix_G_lighting(grad_x, grad_y, p_quad, polygon)

G = zeros(polygon.size);

for i = 1 : polygon.size

    integrand = grad_x(:,:,:,i).^2 + grad_y(:,:,:,i).^2;
    G(i,i) = quadrature_2D(p_quad, integrand,"Evaluated");

    for j = i+1 : polygon.size

        integrand = grad_x(:,:,:,i).*grad_x(:,:,:,j) + grad_y(:,:,:,i).*grad_y(:,:,:,j);
        G(i,j) = quadrature_2D(p_quad, integrand,"Evaluated");
        G(j,i) = G(i,j);

    end

end