function [D] = vem_matrix_D(base, polynomial, polygon, p_quad, base_val) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: vem_matrix_D
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function computes the matrix that contais the DOFs of the basis polynomial function
%
% Input
% ===== 
% base       : Basis function
% polynomial : Information on the polynomial used to compute G
% polygon    : Information on the polygon
%
%
% Output
% ======
% D : The matrix D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

D = zeros(polygon.size, polynomial.dim);

nedges = polygon.nedges;
x_vert = polygon.vertex(1:end-1,1);
y_vert = polygon.vertex(1:end-1,2);
k      = polynomial.k;

for j = 1:polynomial.dim

    D(1:nedges,j) = base(x_vert,y_vert,polynomial.deg(j,:));

end

if (polynomial.k >= 2)
    for j = 1:nedges
    
        %int_x = linspace(polygon.vertex(j,1), polygon.vertex(j+1,1),polynomial.k+1);
        %int_y = linspace(polygon.vertex(j,2), polygon.vertex(j+1,2),polynomial.k+1);
        int_x = polynomial.lobatto * (polygon.vertex(j+1,1)- polygon.vertex(j,1)) ...
              + polygon.vertex(j,1);
        int_y = polynomial.lobatto * (polygon.vertex(j+1,2)- polygon.vertex(j,2)) ...
              + polygon.vertex(j,2);
        int_x = int_x(2:end-1);
        int_y = int_y(2:end-1);

        for d = 1:polynomial.dim

            D(nedges + 1+(j-1)*(k-1): nedges + j*(k-1),d) = base(int_x,int_y,polynomial.deg(d,:));

        end

    end

    for i = 1:polynomial.int

        for j = 1:polynomial.dim
    
            integrand = base_val(:,:,:,i) .* base_val(:,:,:,j);
            D(polynomial.k*polygon.nedges + i,j) = quadrature_2D(p_quad,integrand,"Evaluated") ...
                                                 ./ polygon.area;

        end

    end



end
    
end