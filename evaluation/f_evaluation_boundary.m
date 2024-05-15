function [f_val] = f_evalutation_boundary(f, b_quad)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: f_evalutation_boundary
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function evaluates f on the quadrature nodes in the interior of a polygon
%
% Input
% ===== 
% base   : The function
% p_quad : The quadrature formula
%
% Output
% ======
% f_val : The pointwise values
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% May  07, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

f_val = zeros(size(b_quad.x));                                                                       %Memory allocation
f_val = [f_val(1,:); f_val];

f_val(2:end,:) = f(b_quad.x,b_quad.y);                                   %Evalutation
f_val(1,:)     = f_val(end,:);

end