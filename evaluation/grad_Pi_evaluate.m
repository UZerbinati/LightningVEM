function [Pi_val] = grad_Pi_evaluate(Pi, grad)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: Pi_evaluate
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function evaluates the gradient of the basis function on the quadrature nodes
%
% Input
% ===== 
% Pi       : Coefficients of the projection
% base     : Values of the basis functions
%
% Output
% ======
% Pi_val : Pointwise values of the projection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Pi_val = zeros(size(grad,1:3));                                                                      %Memory allocation
Pi_val(:,:,:, 2*size(Pi,2)) = 0;

dim = size(Pi,2);                                                                                    %Virtual functions
pol = size(Pi,1);                                                                                    %Polygons

for n = 1:size(Pi,2)

    for d = 1:size(Pi,1)
        Pi_val(:,:,:,n) = Pi_val(:,:,:,n) + grad(:,:,:,d) * Pi(d,n);
        Pi_val(:,:,:,n + dim) = Pi_val(:,:,:,n + dim) + grad(:,:,:,d + pol) * Pi(d,n);

    end 

end

end