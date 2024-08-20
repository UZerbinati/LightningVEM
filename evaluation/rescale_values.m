function [val_int, valx_int, valy_int] = rescale_values(val_int, valx_int, valy_int, Coeff, kn)

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

aux   = zeros(size(val_int));
aux_x = zeros(size(val_int));
aux_y = zeros(size(val_int));

for n = 1:size(Coeff,2)

    aux(:,:,:,n)   = sum(val_int(:,:,:,1:kn).*reshape(Coeff(:,n), 1, 1, 1, []),4);
    aux_x(:,:,:,n) = sum(valx_int(:,:,:,1:kn).*reshape(Coeff(:,n), 1, 1, 1, []),4);
    aux_y(:,:,:,n) = sum(valy_int(:,:,:,1:kn).*reshape(Coeff(:,n), 1, 1, 1, []),4);

end

val_int(:,:,:,1:size(Coeff,2))  = aux(:,:,:,1:size(Coeff,2));
valx_int(:,:,:,1:size(Coeff,2)) = aux_x(:,:,:,1:size(Coeff,2));
valy_int(:,:,:,1:size(Coeff,2)) = aux_y(:,:,:,1:size(Coeff,2));

end