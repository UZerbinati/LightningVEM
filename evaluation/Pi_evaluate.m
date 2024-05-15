function [Pi_val] = Pi_evaluate(Pi, base)

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

Pi_val = zeros(size(base,1:3));                                                                      %Memory allocation
Pi_val(:,:,:,size(Pi,2)) = 0;

for n = 1:size(Pi,2)

%     for d = 1:size(base,4)
% 
%         Pi_val(:,:,:,n) = Pi_val(:,:,:,n) + base(:,:,:,d) * Pi(d,n);
% 
%     end 

Pi_val(:,:,:,n) = sum(base.*reshape(Pi(:,n), 1, 1, 1, []),4);
end

end