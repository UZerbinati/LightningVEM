function [int] = quadrature_2D(p_quad, f, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: quadrature_2D 
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function adds the information relative to the edges to the struct domainMesh
%
% Input
% =====
% p_quad  : The quadrature formula
% f       : The function
% options : This variable is set to the default value "NonEvaluated", it means that f is an handle 
%           function, it is possible to set this variable to "Evaluated" if f contains the 
%           pointwise values of f
%
% Output
% ======
% int : The value of the integral
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mai  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

arguments
    p_quad 
    f
    options = "NonEvaluated"

end

if (options == "Evaluated")
    int = sum(p_quad.omega .* f,'all');
    
elseif (options == "NonEvaluated")
    int = 0;
    for i = 1:size(p_quad.omega,1)

        if (p_quad.bool(i) ~= 0) 
            for j = 1:size(p_quad.omega,2)
                for k=1:size(p_quad.omega,3)
                    int = int + p_quad.omega(i,j,k)*f(p_quad.xi(i,j,k),p_quad.eta(i,j));
                end
            end
        end
    end

else
    error("Fai schifo! Punto balzo");
end

end