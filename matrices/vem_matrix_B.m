function [B] = vem_matrix_B(grad_val_bound, polynomial, polygon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_matrix_B
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function computes the matrix that contais the integrals of the scalar products between the 
% gradients of two polynomial basis functions. The first row is replaced by the integral mean on 
% the boundary
%
% Input
% ===== 
% grad_val_int   : Values of the gradients of the basis functions on the boundary
% polynomial     : Information on the polynomial used to compute G
% polygon        : Information on the polygon
%
% Output
% ======
% B : The matrix B
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

B       = zeros(polynomial.dim, polygon.size);

edges   = [polygon.edges(end); polygon.edges];
nedges  = polygon.nedges;
vnormal = [polygon.vnormal(end,:); polygon.vnormal];

weights = polynomial.weights;
val     = polynomial.val;
k       = polynomial.k;

%% PRIMA RIGA
B(1,1:nedges) = (weights * val(:,1)).*(edges(1:end-1) + edges(2:end)) / polygon.perimeter;           %DOFs at the vertices

if (polynomial.k >= 2)                                                                               %2 DOFs
    
    for j = 1:nedges
                    
        B(1,nedges + 1+(j-1)*(k-1): nedges + j*(k-1)) = (weights * val(:,2:end)).* edges(j+1) ...
                                                      / polygon.perimeter;

    end

end

%% RESTO DELLE RIGHE

for j = 1:nedges                                                                                    

    for i = 2:polynomial.dim    

        index = nedges + 1+(j-1)*(k-1): nedges + j*(k-1);
        
        B(i,j) = weights          * ( polynomial.val(:,1) ... 
               .* (vnormal(j,:)   * grad_val_bound([j, j+nedges+1],end:-1:1,i))' ) * edges(j);
        B(i,j) = B(i,j) + weights * (polynomial.val(:,1) ...
               .* (vnormal(j+1,:) * grad_val_bound([j+1, j+2+nedges],:,i))'      ) * edges(j+1);
        B(i,index) = weights * (polynomial.val(:,2:end) ...
               .* (vnormal(j+1,:) * grad_val_bound([j+1, j+2+nedges],:,i))'      ) * edges(j+1);

    end

end    

for j = 1:polynomial.dim                                                                             %DOF 3 INTERNO

    deg    = polynomial.deg(j,:);
    deg_xx = -1;
    deg_yy = -1;

    if ( deg(1) >= 2 )

        deg_xx = [deg(1)-2, deg(2)];

    end

    if ( deg(2) >= 2 )

        deg_yy = [deg(1), deg(2)-2];

    end

    for d = 1:polynomial.int

        if (polynomial.deg(d,:) == deg_xx)

            B(j,nedges*k + d) = B(j,nedges*k + d) - (deg(1)^2 - deg(1))*polygon.area / polygon.diameter / polygon.diameter;

        end

        if (polynomial.deg(d,:) == deg_yy)

            B(j,nedges*k + d) = B(j,nedges*k + d) - (deg(2)^2 - deg(2)) * polygon.area / polygon.diameter / polygon.diameter;

        end

    end
    
end

end