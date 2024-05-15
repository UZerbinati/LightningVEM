function [vem_val] = vem_evalutation_interior(siz, u_vem, p_quad)

% vem_evalutation_interior: this functions evaluates the vem basis function on
%                           the quadrature nodes.
%
% Input parameters:
%    siz: number of basis functions;
%  u_vem: the basis functions;
% p_quad: quadrature points.
%
% Output parameters:
% vem_val: pointwise values.

vem_val = zeros(size(p_quad.xi));

vem_val(:,:,:,siz) = 0;

for j = 1 : siz

    vem_val(:,:,:,j) = u_vem{j}(p_quad.xi + 1i.*p_quad.eta);

end

end