function [errL2, errH1, values] = compute_errors(domainMesh, polynomial, U, u, deru_x, deru_y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: compute_errors
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function returns the errors in the L2 norm and the H1 seminorm
%
% Input
% =====
% domainMesh : The struct that contains the information on the mesh that we are using
% polynomial : The struct that contains the information on the polynomials that we are using
% U          : The numerical solution vector
% u          : The analytic solution
% deru_x     : x-derivative of the solution
% deru_y     : y-derivative of the solution
%
% Output
% ======
% errL2 : L2 norm of the error
% errH1 : H1 seminorm of the error
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% May  07, 2022: first realease (by M. Trezzi)
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    errL2 = 0;
    errH1 = 0;
    values = zeros(domainMesh.npolygon,1);

    for i=1:domainMesh.npolygon

        polygon = domainMesh.polygon{i};                                                             %Information on the i-th polygon

        [p_quad] = polygon_quadrature(domainMesh.polygon{i},polynomial.n+1);                         %Quadrature formula 

        base = @(x,y,k) ((x - polygon.centroid(1)) ./ polygon.diameter).^k(:,1) ...
             .* ((y - polygon.centroid(2)) ./ polygon.diameter).^k(:,2);
        grad = @(x,y,k) grad_fun(x,y,k,polygon);

        base_val   = base_evalutation_interior(p_quad, polynomial, domainMesh.polygon{i});           %Evalutations of the basis function
        grad_val   = grad_evalutation_interior(base_val, p_quad, polynomial, domainMesh.polygon{i});

        u_val      = f_evalutation_interior(u,      p_quad);                                         %Evalutations of the analytic solution
        deru_x_val = f_evalutation_interior(deru_x, p_quad);
        deru_y_val = f_evalutation_interior(deru_y, p_quad);

        Pistar = polygon.G \ (polygon.B * U(polygon.local_dofs));                                    %Coefficients of the Pistar projection

        derPistar_val_x = Pi_evaluate(Pistar, grad_val(:,:,:,1:polynomial.dim));
        derPistar_val_y = Pi_evaluate(Pistar, grad_val(:,:,:,polynomial.dim+1:end));

        P0 = polygon.H \ (polygon.C * U(polygon.local_dofs));                                        %Coefficients of the Pistar projection

        P0_val = Pi_evaluate(P0, base_val);                                                          %Evalutations of the Pistar projection
        
        diff_sq      = (u_val - P0_val).^2;                                                      %Integrands
        diff_sq_grad = (deru_x_val - derPistar_val_x).^2 ...
                     + (deru_y_val - derPistar_val_y).^2;

        errL2 = errL2 + quadrature_2D(p_quad,diff_sq,     "Evaluated");                              %Update errors                                    
        errH1 = errH1 + quadrature_2D(p_quad,diff_sq_grad,"Evaluated");

        values(i) = quadrature_2D(p_quad,  P0_val.^2, "Evaluated");

    end

    errL2 = sqrt(errL2);
    errH1 = sqrt(errH1);

end