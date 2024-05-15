function [errL2, errH1] = compute_errors_lighting(domainMesh, U, u, deru_x, deru_y, k)

% compute_errors_lighting: compute the errors in the L^2 norm and H^2 seminorm.
%
% Input parameters:
% domainMesh: the structure that contains the information on the mesh;
%          U: vector that contains the numerical solution;
%          u: continuous solution;
%     deru_x: x-derivative of u;
%     deru_y: y-derivative of u;
%          k: degree of the polynomials.
%
% Output parameter:
% errL2: error in the L2 norm;
% errH1: error in the H1 seminorm.

polygons = domainMesh.polygon;

errL2vec = zeros(domainMesh.npolygon,1);
errH1vec = zeros(domainMesh.npolygon,1);

U_vec    = parallel.pool.Constant(U);

parfor i=1:domainMesh.npolygon

    polygon = polygons{i};                                                                           %Information on the i-th polygon

    [p_quad] = polygon_quadrature(polygon,k+3);                                                      %Quadrature formula 

    u_val      = f_evaluation_interior(u,      p_quad);                                              %evaluations of the analytic solution
    deru_x_val = f_evaluation_interior(deru_x, p_quad);
    deru_y_val = f_evaluation_interior(deru_y, p_quad);

    u_sol_val       = zeros(size(p_quad.xi));
    val_int_xp  = zeros(size(p_quad.xi));
    val_int_xm  = zeros(size(p_quad.xi));
    val_int_yp  = zeros(size(p_quad.xi));
    val_int_ym  = zeros(size(p_quad.xi));
    val_int_2xp = zeros(size(p_quad.xi));
    val_int_2xm = zeros(size(p_quad.xi));
    val_int_2yp = zeros(size(p_quad.xi));
    val_int_2ym = zeros(size(p_quad.xi));

    h = 1e-6*polygon.diameter;                                                                       %Size of the finite difference

    U_val = U_vec.Value(polygon.local_dofs);

    % Approximate the derivative of the analytic solution with a second order finite difference
    for j = 1:polygon.size
    
        u_sol_val = u_sol_val + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*p_quad.eta);

        val_int_xp = val_int_xp + U_val(j) * polygon.basis{j}(p_quad.xi + h + 1i.*p_quad.eta);

        val_int_xm = val_int_xm + U_val(j) * polygon.basis{j}(p_quad.xi - h + 1i.*p_quad.eta);
        
        val_int_yp = val_int_yp + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*(p_quad.eta + h));

        val_int_ym = val_int_ym + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*(p_quad.eta - h));

        %+2h

        val_int_2xp = val_int_2xp + U_val(j) * polygon.basis{j}(p_quad.xi + 2*h + 1i.*p_quad.eta);
        val_int_2xm = val_int_2xm + U_val(j) * polygon.basis{j}(p_quad.xi - 2*h + 1i.*p_quad.eta);
        val_int_2yp = val_int_2yp + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*(p_quad.eta + 2*h));
        val_int_2ym = val_int_2ym + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*(p_quad.eta - 2*h));
        
    end

    diff_sq     = (u_sol_val - u_val).^2;                                                            %Difference u - u_h
    errL2vec(i) = quadrature_2D(p_quad,diff_sq,"Evaluated");

    grad_val_int_x = (val_int_xp - val_int_xm) / (2*h);
    grad_val_int_y = (val_int_yp - val_int_ym) / (2*h);
    
    diff_sq     = (deru_x_val - grad_val_int_x).^2 + (grad_val_int_y - deru_y_val).^2;               %grad(u-u_h)
    errH1vec(i) = quadrature_2D(p_quad,diff_sq,"Evaluated");

end

errL2 = sqrt(sum(errL2vec));
errH1 = sqrt(sum(errH1vec));

end
