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

    [p_quad] = polygon_quadrature(polygon,k+12);                                                      %Quadrature formula 

    u_val      = f_evaluation_interior(u,      p_quad);                                              %evaluations of the analytic solution
    deru_x_val = f_evaluation_interior(deru_x, p_quad);
    deru_y_val = f_evaluation_interior(deru_y, p_quad);

    u_sol_val  = zeros(size(u_val));
    u_sol_valx = zeros(size(u_val));
    u_sol_valy = zeros(size(u_val));

    val_int  = vem_evaluation_interior(polygon.size, polygon.basis,  p_quad);                                      %Evaluates the basis function in the quadrature nodes
    valx_int = vem_evaluation_interior(polygon.size, polygon.basisX, p_quad);
    valy_int = vem_evaluation_interior(polygon.size, polygon.basisY, p_quad);

    [val_int, valx_int, valy_int] = rescale_values(val_int, valx_int, valy_int, polygon.coeff, size(polygon.coeff,1));

    if k == 3
        range = k*polygon.nedges + (1:3);

        [val_int(:,:,:,range), valx_int(:,:,:,range), valy_int(:,:,:,range)] = rescale_values(val_int(:,:,:,range), valx_int(:,:,:,range), valy_int(:,:,:,range), polygon.coeff2, 3);
    end

    U_val = U_vec.Value(polygon.local_dofs);

    % Approximate the derivative of the analytic solution with a second order finite difference
    for j = 1:polygon.size
    
%         u_sol_val  = u_sol_val  + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*p_quad.eta);
%         u_sol_valx = u_sol_valx + U_val(j) * polygon.basisX{j}(p_quad.xi + 1i.*p_quad.eta);
%         u_sol_valy = u_sol_valy + U_val(j) * polygon.basisY{j}(p_quad.xi + 1i.*p_quad.eta);

        u_sol_val  = u_sol_val  + U_val(j) * val_int(:,:,:,j);
        u_sol_valx = u_sol_valx + U_val(j) * valx_int(:,:,:,j);
        u_sol_valy = u_sol_valy + U_val(j) * valy_int(:,:,:,j);
        
    end

    diff_sq     = (u_sol_val - u_val).^2;                                                            %Difference u - u_h
    errL2vec(i) = max(quadrature_2D(p_quad, diff_sq, "Evaluated"),0);
    
    %diff_sq     = (deru_x_val - grad_val_int_x).^2 + (grad_val_int_y - deru_y_val).^2;               %grad(u-u_h)
    diff_sq     = (deru_x_val - u_sol_valx).^2 + (deru_y_val - u_sol_valy).^2;
    errH1vec(i) = max(quadrature_2D(p_quad,diff_sq,"Evaluated"),0);

end

errL2 = sqrt(sum(errL2vec));
errH1 = sqrt(sum(errH1vec));

end
