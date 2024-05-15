function [base_val] = base_evalutation_boundary(b_quad, polynomial, polygon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: base_evalutation_boundary
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function evaluates the basis function on the quadrature nodes on the boundary of a polygon
%
% Input
% ===== 
% b_quad     : The quadrature formula
% polynomial : The struct that contains the information on the polynomials that we are using
%
% Output
% ======
% base_val : The pointwise values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

dims   = size(b_quad.x);
dims(1) = dims(1) + 1;
base_val = ones(dims);                                                                    %Memory allocation
%base_val = [base_val(1,:); base_val];

base_val(:,:,polynomial.dim) = 1;

diffx = (b_quad.x - polygon.centroid(1)) ./ polygon.diameter;
diffy = (b_quad.y - polygon.centroid(2)) ./ polygon.diameter;

for d = 2:polynomial.dim                                                                             %Evalutation

    %base_val(2:end,:,d) = base(b_quad.x,b_quad.y,polynomial.deg(d,:));
    if polynomial.deg(d,1) ~= 0
        base_val(2:end,:,d) = base_val(2:end,:,d-sum(polynomial.deg(d,:))) .* diffx;
    else
        base_val(2:end,:,d) = base_val(2:end,:,d-sum(polynomial.deg(d,:))-1) .* diffy;
    end
    
end

base_val(1,:,:) = base_val(end,:,:);

end