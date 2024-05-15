function [G] = vem_matrix_G(base_val_bound, grad_val_int, p_quad, polynomial, polygon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_matrix_G
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function computes the matrix that contais the integrals of the scalar products between the 
% gradients of a polynomial and a virtual element function. The first row
% is the integral mean over the boundary
%
% Input
% ===== 
% base_val_bound : Values of the basis functions on the boundary
% grad_val_int   : Values of the gradients of the basis functions in the interior
% p_quad         : Quadrature formula
% polynomial     : Information on the polynomial used to compute G
% polygon        : Information on the polygon
%
%
% Output
% ======
% G : The matrix G
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

dim = polynomial.dim;                                                                                %Memory allocation
G   = zeros(dim, dim);

for i=2:dim

    integrand = grad_val_int(:,:,:,i).^2 + grad_val_int(:,:,:,i+dim).^2;                             %Diagonal terms
    
    G(i,i) = quadrature_2D(p_quad, integrand, "Evaluated");

    for j=i+1:dim

       integrand = grad_val_int(:,:,:,i)     .* grad_val_int(:,:,:,j) ...                            %Upper triangluar terms
                 + grad_val_int(:,:,:,i+dim) .* grad_val_int(:,:,:,j+dim);

       G(i,j) = quadrature_2D(p_quad, integrand, "Evaluated");
       G(j,i) = G(i,j);

    end

end

G(1,:) = projection_P0(base_val_bound, polynomial, polygon);                                         %Computes the first row

end