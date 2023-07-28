function [f_val] = f_evalutation_interior(f, p_quad)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: f_evalutation_interior
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

f_val = f(p_quad.xi,p_quad.eta);

end