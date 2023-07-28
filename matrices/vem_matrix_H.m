function [H] = vem_matrix_H(base_val_int, p_quad, polynomial)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_matrix_H
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function computes the matrix that contais the integrals of the products of two polynomial 
% basis function
%
% Input
% ===== 
% base_val_int : Values of the basis functions in the interior
% p_quad       : Quadrature formula
% polynomial   : Information on the polynomial used to compute G
%
% Output
% ======
% H : The matrix H
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% May.  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

H = zeros(polynomial.dim, polynomial.dim);

for i=1:polynomial.dim

    for j=i:polynomial.dim
        
        integrand = base_val_int(:,:,:,i) .* base_val_int(:,:,:,j);

        H(i,j) = quadrature_2D(p_quad, integrand, "Evaluated");
        H(j,i) = H(i,j);
        
    end

end

end