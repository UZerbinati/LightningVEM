function [C] = vem_matrix_C(base_val_int, Pi, p_quad, polynomial, polygon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_matrix_C
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function computes the matrix that contains the integrals of a virtual function times a 
% polynomial
%
% Input
% ===== 
% base_val_int : Values of the basis functions in the interior
% Pi           : Coefficients of the Pistar projection
% p_quad       : Quadrature formula
% polynomial   : Information on the polynomial used to compute G
% polygon      : Information on the polygon
%
% Output
% ======
% C : The matrix C
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

C = zeros(polynomial.dim, polygon.size);

% for i = 1:polynomial.dim
%     for j = 1:polygon.size
%         
%             integrand = base_val_int(:,:,:,i) .* Pi(:,:,:,j);
% 
%             C(i,j) = quadrature_2D(p_quad, integrand, "Evaluated");
%            
%     end
% end

index = (polygon.nedges * polynomial.k + 1:1:polygon.size);
C(1:polynomial.int ,index) = diag(ones(polynomial.int,1)).*polygon.area;

for i = polynomial.int + 1:polynomial.dim
    for j = 1:polygon.size
        
            integrand = base_val_int(:,:,:,i) .* Pi(:,:,:,j);

            C(i,j) = quadrature_2D(p_quad, integrand, "Evaluated");
           
    end
end

end


