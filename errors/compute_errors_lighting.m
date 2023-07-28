function [errL2, errH1] = compute_errors_lighting(domainMesh, U, u, deru_x, deru_y, k, hh)

polygons = domainMesh.polygon;

errL2vec = zeros(domainMesh.npolygon,1);
errH1vec = zeros(domainMesh.npolygon,1);

U_vec    = parallel.pool.Constant(U);

parfor i=1:domainMesh.npolygon

    polygon = polygons{i};                                                                           %Information on the i-th polygon

    [p_quad] = polygon_quadrature(polygon,k+3);                                                      %Quadrature formula 

    u_val      = f_evalutation_interior(u,      p_quad);                                             %Evalutations of the analytic solution
    deru_x_val = f_evalutation_interior(deru_x, p_quad);
    deru_y_val = f_evalutation_interior(deru_y, p_quad);

    u_sol_val       = zeros(size(p_quad.xi));
    vem_val_int_xp  = zeros(size(p_quad.xi));
    vem_val_int_xm  = zeros(size(p_quad.xi));
    vem_val_int_yp  = zeros(size(p_quad.xi));
    vem_val_int_ym  = zeros(size(p_quad.xi));
    vem_val_int_2xp = zeros(size(p_quad.xi));
    vem_val_int_2xm = zeros(size(p_quad.xi));
    vem_val_int_2yp = zeros(size(p_quad.xi));
    vem_val_int_2ym = zeros(size(p_quad.xi));

    h = hh*polygon.diameter;

    U_val = U_vec.Value(polygon.local_dofs);

    for j = 1:polygon.size
    
        u_sol_val = u_sol_val + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*p_quad.eta);

        vem_val_int_xp = vem_val_int_xp + U_val(j) * polygon.basis{j}(p_quad.xi + h + 1i.*p_quad.eta);

        vem_val_int_xm = vem_val_int_xm + U_val(j) * polygon.basis{j}(p_quad.xi - h + 1i.*p_quad.eta);
        
        vem_val_int_yp = vem_val_int_yp + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*(p_quad.eta + h));

        vem_val_int_ym = vem_val_int_ym + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*(p_quad.eta - h));

        %+2h

        vem_val_int_2xp = vem_val_int_2xp + U_val(j) * polygon.basis{j}(p_quad.xi + 2*h + 1i.*p_quad.eta);
        vem_val_int_2xm = vem_val_int_2xm + U_val(j) * polygon.basis{j}(p_quad.xi - 2*h + 1i.*p_quad.eta);
        vem_val_int_2yp = vem_val_int_2yp + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*(p_quad.eta + 2*h));
        vem_val_int_2ym = vem_val_int_2ym + U_val(j) * polygon.basis{j}(p_quad.xi + 1i.*(p_quad.eta - 2*h));
        
    end

    diff_sq = (u_sol_val - u_val).^2;
    errL2vec(i) = quadrature_2D(p_quad,diff_sq,"Evaluated");

    grad_val_int_x = (vem_val_int_xp - vem_val_int_xm) / (2*h);
    grad_val_int_y = (vem_val_int_yp - vem_val_int_ym) / (2*h);
    
    diff_sq = (deru_x_val - grad_val_int_x).^2 + (grad_val_int_y - deru_y_val).^2;
    errH1vec(i) = quadrature_2D(p_quad,diff_sq,"Evaluated");

end

errL2 = sqrt(sum(errL2vec));
errH1 = sqrt(sum(errH1vec));

end
