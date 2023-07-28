function [vem_val] = vem_evalutation_interior(siz, u_vem, p_quad)

vem_val = zeros(size(p_quad.xi));

vem_val(:,:,:,siz) = 0;

for j = 1 : siz

    vem_val(:,:,:,j) = u_vem{j}(p_quad.xi + 1i.*p_quad.eta);

end

end